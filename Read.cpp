#include "Read.h"


void Read::kmersOf(void) {
    for (int i = 0; i < sequence.size() - kmerLength; ++i) {
	string kmer = sequence.substr(i, kmerLength);
	(*kmersToReads)[kmer].push_back(this);
	map<string, vector<Read*> >::iterator k = kmersToReads->find(kmer);
	kmers.push_back(k);
    }
}


void Read::buildMutualKmerHits(void) {
    for (vector<map<string, vector<Read*> >::iterator>::iterator i = kmers.begin(); i != kmers.end(); ++i) {
	vector<Read*>& reads = (*i)->second;
	for (vector<Read*>::iterator r = reads.begin(); r != reads.end(); ++r) {
	    ++mutualKmerHits[*r];
	}
    }
    for (map<Read*, int>::iterator m = mutualKmerHits.begin(); m != mutualKmerHits.end(); ++m) {
	mutualKmerHitCounts[m->second].push_back(m->first);
    }
}

void Read::alignWithReads(vector<Read*>& reads, bool useSSW) {
    for (vector<Read*>::iterator r = reads.begin(); r != reads.end(); ++r) {
	Read* read = *r;
	if (read == this) {
	    continue;
	}
	if (useSSW) {
	    aligner.Align(read->sequence.c_str(), sequence.c_str(), sequence.size(), filter, &alignment);
	    alignments[read] = alignment;
	    swScores[alignment.sw_score].push_back(read);
	    int alignmentLength = alignment.query_end - alignment.query_begin;
	    float normalizedSwScore = (float) alignment.sw_score / (float) alignmentLength;
	    swScoresPerBase[normalizedSwScore].push_back(read);
	} else {	
	    CSmithWatermanGotoh sw(2, 2, 3, 1); // testing only, uses the same values as SSW
	    unsigned int referencePos;
	    string cigar;
	    float bestScore = 0;
	    sw.Align(referencePos, cigar, sequence, read->sequence);
	    alignment.ref_begin = referencePos;
	    alignment.sw_score = sw.BestScore;
	    alignment.cigar_string = cigar;
	    alignments[read] = alignment;
	    swScores[alignment.sw_score].push_back(read);
	    int alignmentLength = sequence.size();
	    float normalizedSwScore = (float) alignment.sw_score / (float) alignmentLength;
	    swScoresPerBase[normalizedSwScore].push_back(read);
	}
    }
}

void Read::alignWithCandidates(int minKmerMatchCount, bool useSSW) {
    for (map<int, vector<Read*> >::iterator rs = mutualKmerHitCounts.begin(); rs != mutualKmerHitCounts.end(); ++rs) {
	if (rs->first < minKmerMatchCount) {
	    continue;
	}
	alignWithReads(rs->second, useSSW);
    }
}

ostream& operator<<(ostream& o, Read& read) {
    int begin = 0;
    for (map<Read*, StripedSmithWaterman::Alignment>::iterator m = read.alignments.begin();
	 m != read.alignments.end(); ++m) {
	StripedSmithWaterman::Alignment& alignment = m->second;
	//Read& oread = *m->first;
	if (begin < alignment.query_begin) {
	    begin = alignment.query_begin;
	}
    }

//    for (map<Read*, StripedSmithWaterman::Alignment>::iterator m = read.alignments.begin();
//	 m != read.alignments.end(); ++m) {

    int startoffset = 0;
    for (map<int, vector<Read*> >::iterator p = read.swScores.begin(); p != read.swScores.end(); ++p) {
	for (vector<Read*>::iterator r = p->second.begin(); r != p->second.end(); ++r) {
	    StripedSmithWaterman::Alignment& alignment = read.alignments[*r];
	    vector<pair<int, char> > scigar = splitCigar(alignment.cigar_string);
	    int sclips = (scigar.front().second == 'S' ? scigar.front().first : 0 );
	    int start = alignment.ref_begin - sclips;
	    if (start < startoffset) {
		startoffset = start;
	    }
	}
    }

    int padding = abs(startoffset) + read.name.size() + 40;
    o << read.name << string(padding - read.name.size(), ' ') << read.sequence << endl;

    for (map<float, vector<Read*> >::iterator p = read.swScoresPerBase.begin(); p != read.swScoresPerBase.end(); ++p) {
	float scorePerBase = p->first;
	for (vector<Read*>::iterator r = p->second.begin(); r != p->second.end(); ++r) {
	    StripedSmithWaterman::Alignment& alignment = read.alignments[*r];
	    Read& oread = **r;

	    vector<pair<int, char> > scigar = splitCigar(alignment.cigar_string);
	    int sclipsStart = (scigar.front().second == 'S' ? scigar.front().first : 0 );
	    int sclipsEnd = (scigar.back().second == 'S' ? scigar.back().first : 0 );
	    int start = alignment.ref_begin - sclipsStart;

	    // deal with soft clips
	    string adjustedSequence = oread.sequence;
	    //transform(adjustedSequence.begin(), adjustedSequence.begin() + sclipsStart, adjustedSequence.begin(), ::tolower);
	    //transform(adjustedSequence.rbegin(), adjustedSequence.rbegin() + sclipsEnd, adjustedSequence.rbegin(), ::tolower);

	    // deal with indels in the same fashion
	    int rpos = 0;
	    for (vector<pair<int, char> >::iterator s = scigar.begin(); s != scigar.end(); ++s) {
		int l = s->first;
		char c = s->second;
		switch (c) {
		case 'S':
		case 'I':
		    transform(adjustedSequence.begin() + rpos,
			      adjustedSequence.begin() + rpos + l,
			      adjustedSequence.begin() + rpos,
			      ::tolower);
		    rpos += l;
		    break;
		case 'M':
		    rpos += l;
		    break;
		case 'D':
		    break;
		default:
		    break;
		}
	    }

	    stringstream ss;
	    ss << scorePerBase << ":" << alignment.sw_score << ":" << alignment.ref_begin
	       << ":" << alignment.cigar_string << " " << oread.name;
	    string nameetc = ss.str();

	    o << nameetc << string(max(1, (int) (padding + start - nameetc.size())), ' ') << adjustedSequence << endl;
	}
    }
    o << endl;
    return o;
}

void Read::dotFormat(ostream& out, float minPerBaseSW, int swScoreCut) {
    for (map<float, vector<Read*> >::iterator p = swScoresPerBase.begin(); p != swScoresPerBase.end(); ++p) {
	if (p->first < minPerBaseSW) continue;
	for (vector<Read*>::iterator r = p->second.begin(); r != p->second.end(); ++r) {
	    StripedSmithWaterman::Alignment& alignment = alignments[*r];
	    Read& read = **r;
	    if (alignment.sw_score > swScoreCut) {
		if (alignment.ref_begin > 0) {
		    out << "\"" << name << "\"" << " -> " << "\"" << read.name << "\""
			<< " [ weight=\"" << p->first << "\",label=\"" << p->first << "\"];" << endl;
		    //<< " [ weight=\"" << alignment.sw_score << "\",label=\"" << alignment.sw_score << "\"];" << endl;
		} else {
		    out << "\"" << read.name << "\"" << " -> " << "\"" << name << "\""
			<< " [ weight=\"" << p->first << "\",label=\"" << p->first << "\"];" << endl;
		}
	    }
	}
    }
}

ostream& operator<<(ostream& out, vector<pair<int, char> >& cigar) {
    for (vector<pair<int, char> >::iterator c = cigar.begin(); c != cigar.end(); ++c) {
	out << c->first << c->second;
    }
    return out;
}

vector<pair<int, char> > splitCigar(const string& cigarStr) {
    vector<pair<int, char> > cigar;
    string number;
    char type = '\0';
    // strings go [Number][Type] ...
    for (string::const_iterator s = cigarStr.begin(); s != cigarStr.end(); ++s) {
        char c = *s;
        if (isdigit(c)) {
            if (type == '\0') {
                number += c;
            } else {
                // signal for next token, push back the last pair, clean up
                cigar.push_back(make_pair(atoi(number.c_str()), type));
                number.clear();
		type = '\0';
                number += c;
            }
        } else {
            type = c;
        }
    }
    if (!number.empty() && type != '\0') {
        cigar.push_back(make_pair(atoi(number.c_str()), type));
    }
    return cigar;
}
