

// plan
// 
// we need a number of items:
// - phased vcf to haplotypes (in vcf) ==> fasta
// - haplotypes to reads (wgsim)  (in fasta) ==> fastq
// - align rp. to ref
// - read assembly via alignment (reads)
//   * map<read-id, map<sswscore, alignment*> > -- provides traversal of graph
//   * map<read-id, map<kmer-mutual-hit-count, alignment*> > -- determines reads to align
//   * map<kmer, vector<read-id> > -- use to generate kmer-mutual-hit-count
//   * once structures have been generated, run through each kmer-mutual-hit-count aligning against reads when
//     @ they haven't been aligned
//     @ the alignment meets some threshold of quality (low, but enough that the read must align partly)
//     @ record alignments meeting some threshold of quality
//   * walk back across table, building connection graph between reads

#include <iostream>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <signal.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <map>
#include <vector>

#include "fastahack/Fasta.h"

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAlignment.h"

#include "ssw_cpp.h"
#include "Read.h"

using namespace std;
using namespace BamTools;


short qualityChar2ShortInt(char c) {
    return static_cast<short>(c) - 33;
}



void countMismatchesAndGaps(
    BamAlignment& alignment,
    vector<CigarOp>& cigarData,
    string referenceSequence,
    int& mismatches,
    int& gaps,
    int& gapslen,
    int& softclips,
    int& mismatchQsum,
    int& softclipQsum
    ) {

    int sp = 0;
    int rp = 0;
    for (vector<CigarOp>::const_iterator c = cigarData.begin();
	 c != cigarData.end(); ++c) {
        int l = c->Length;
        char t = c->Type;
        if (t == 'M') { // match or mismatch
            for (int i = 0; i < l; ++i) {
                if (alignment.QueryBases.at(rp) != referenceSequence.at(sp)) {
                    ++mismatches;
		    mismatchQsum += qualityChar2ShortInt(alignment.Qualities.at(rp));
		}
                ++sp;
                ++rp;
            }
        } else if (t == 'D') { // deletion
            ++gaps;
	    gapslen += l;
            sp += l;  // update reference sequence position
        } else if (t == 'I') { // insertion
	    ++gaps;
	    gapslen += l;
	    rp += l;  // update read position
	} else if (t == 'S') { // soft clip, clipped sequence present in the read not matching the reference
	    softclips += l;
	    for (int i = 0; i < l; ++i) {
		softclipQsum += qualityChar2ShortInt(alignment.Qualities.at(rp));
		++rp;
	    }
	} else if (t == 'H') { // hard clip on the read, clipped sequence is not present in the read
	} else if (t == 'N') { // skipped region in the reference not present in read, aka splice
	    sp += l;
	}
    }

}


void printUsage(char** argv) {
    cerr << "usage: [BAM data stream] | " << argv[0] << " [options]" << endl
	 << endl
	 << "Uses alignment to build a string graph between reads in input." << endl
	 << endl
	 << "arguments:" << endl
	 << "    -f --fasta-reference FILE  FASTA reference file to use for realignment (required)" << endl
	 << "    -s --sequence SEQ          Additional sequence to include in mutual alignments." << endl
	 << "    -l --kmer-length N         For pruning possible alignments, use kmer matches of this length." << endl
	 << "    -c --kmer-count N          Only attempt alignments when they share this many kmers matches." << endl
	 << "    -C --stdin                 Read BAM stream from stdin." << endl
	 << "    -Q --soft-clip-qsum N      Trigger assembly if the sum of quality scores of" << endl
	 << "                               soft clipped bases is >= N" << endl
	 << "    -m --match-score N         Set alignment match score. (default 2)" << endl
	 << "    -M --mismatch-penalty N    Set alignment mismatch penalty. (default 2)" << endl
	 << "    -g --gap-open-penalty N    Set alignment gap open penalty. (default 3)" << endl
	 << "    -e --gap-extend-penalty N  Set alignment gap extend penalty. (default 1)" << endl
	 << "    -D --dot-output N          Output a graphviz dot file with per-base SW min of N. (default 1)" << endl
	 << "    -S --dot-output-sw-min N   Only output edges for graphviz dot file if SW score is greater than N. (default 0)" << endl
	 << "    -O --use-standard-sw       Do not use the (very fast) SSW pairwise alignment algorithm." << endl;
}

int main(int argc, char** argv) {

    int c;

    FastaReference reference;
    bool has_ref = false;
    bool debug = false;

    bool stdin = false;
    bool suppress_output = true;

    int kmerLength = 11;
    int kmerCount = 4;

    int matchScore = 2;
    int mismatchPenalty = 2;
    int gapOpenPenalty = 3;
    int gapExtendPenalty = 1;

    bool dotOutput = false;
    float dotPerBaseSWMin = 1;
    int dotSWMin = 0;

    bool useSSW = true;

    int minSoftClipQSum = 20;

    int flanking_window = 50;

    vector<string> sequences;

    if (argc < 2) {
        printUsage(argv);
        exit(1);
    }

    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"debug", no_argument, 0, 'd'},
            {"fasta-reference", required_argument, 0, 'f'},
	    {"kmer-length", required_argument, 0, 'l'},
	    {"kmer-count", required_argument, 0, 'c'},
	    {"stdin", no_argument, 0, 'C'},
	    {"sequence", required_argument, 0, 's'},
	    {"match-score", required_argument, 0, 'm'},
	    {"mismatch-penalty", required_argument, 0, 'M'},
	    {"gap-open-penalty", required_argument, 0, 'g'},
	    {"gap-extend-penalty", required_argument, 0, 'e'},
	    {"dot-output", required_argument, 0, 'D'},
	    {"dot-output-sw-min", required_argument, 0, 'S'},
	    {"soft-clip-qsum", required_argument, 0, 'Q'},
	    {"use-standard-sw", no_argument, 0, 'O'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hdOCf:l:c:s:g:e:m:M:D:S:Q:",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c) {

            case 'f':
                reference.open(optarg); // will exit on open failure
                has_ref = true;
                break;
     
            case 'd':
                debug = true;
                break;

	    case 'l':
		kmerLength = atoi(optarg);
		break;

	    case 'c':
		kmerCount = atoi(optarg);
		break;

	    case 'm':
		matchScore = atoi(optarg);
		break;

	    case 'M':
		mismatchPenalty = atoi(optarg);
		break;

	    case 'g':
		gapOpenPenalty = atoi(optarg);
		break;

	    case 'e':
		gapExtendPenalty = atoi(optarg);
		break;

	    case 'C':
	        stdin = true;
		break;

    	    case 'Q':
		minSoftClipQSum = atoi(optarg);
		break;

	    case 'D':
		dotOutput = true;
		dotPerBaseSWMin = atof(optarg);
		break;

	    case 'S':
		dotSWMin = atoi(optarg);
		break;
	    
	    case 's':
		sequences.push_back(optarg);
		break;

   	    case 'O':
		useSSW = false;
		break;

	    case 'h':
		printUsage(argv);
		exit(1);
		break;

            case '?':
                printUsage(argv);
                exit(1);
                break;
     
	    default:
                abort();
                break;
        }
    }

    if (!has_ref) {
        cerr << "no FASTA reference provided" << endl;
        exit(1);
    }

    BamReader reader;
    if (stdin) {
	if (!reader.Open("stdin")) {
	    cerr << "could not open stdin for reading" << endl;
	    exit(1);
	}
    }

    BamWriter writer;
    if (!suppress_output && !writer.Open("stdout", reader.GetHeaderText(), reader.GetReferenceData())) {
        cerr << "could not open stdout for writing" << endl;
        exit(1);
    }

    // store the names of all the reference sequences in the BAM file
    map<int, string> referenceIDToName;
    vector<RefData> referenceSequences = reader.GetReferenceData();
    int i = 0;
    for (RefVector::iterator r = referenceSequences.begin(); r != referenceSequences.end(); ++r) {
        referenceIDToName[i] = r->RefName;
        ++i;
    }

    BamAlignment alignment;
    //map<long unsigned int, vector<BamAlignment> > alignmentSortQueue;

    map<string, BamAlignment> alignments;
    vector<Read*> reads;
    map<string, vector<Read*> > kmersToReads;

    // for now, operate on whole region
    // next up, tear-down and build-up

    // TODO for assembly, trim repeats from ends of reads (somehow)

    while (reader.GetNextAlignment(alignment)) {
	// get the overlapping reference sequnce to determine mismatches
	if (alignment.IsMapped()) {
	    int endpos = alignment.GetEndPosition();
	    int length = endpos - alignment.Position;  // 0-based half-open interval
	    string ref = reference.getSubSequence(referenceIDToName[alignment.RefID],
						  max(0, alignment.Position - flanking_window),
						  length + 2 * flanking_window);
	    int mismatchesBefore = 0;
	    int gapsBefore = 0;
	    int gapslenBefore = 0;
	    int softclipsBefore = 0;
	    int mismatchQsumBefore = 0;
	    int softclipQsumBefore = 0;
	    countMismatchesAndGaps(alignment,
				   alignment.CigarData,
				   ref.substr(flanking_window, length),
				   mismatchesBefore,
				   gapsBefore,
				   gapslenBefore,
				   softclipsBefore,
				   mismatchQsumBefore,
				   softclipQsumBefore);
	    if (softclipQsumBefore < minSoftClipQSum) {
		continue;
	    }
	}
	alignments[alignment.Name] = alignment;
	BamAlignment* alignmentPtr = &(alignments.find(alignment.Name)->second);
	Read* read = new Read(alignment.QueryBases,
			      &kmersToReads,
			      kmerLength,
			      alignmentPtr,
			      matchScore,
			      mismatchPenalty,
			      gapOpenPenalty,
			      gapExtendPenalty);
	reads.push_back(read);
    }

    BamAlignment emptyBamAlignment;
    emptyBamAlignment.Name = "unknown";

    for (vector<string>::iterator s = sequences.begin(); s != sequences.end(); ++s) {
	Read* read = new Read(*s, &kmersToReads, kmerLength, &emptyBamAlignment);
	reads.push_back(read);
    }

    /*
    for (map<string, vector<Read*> >::iterator k = kmersToReads.begin(); k != kmersToReads.end(); ++k) {
	cout << k->first << endl;
	for (vector<Read*>::iterator r = k->second.begin(); r != k->second.end(); ++r) {
	    cout << "\t" << *r << endl;
	}
    }
    */

    for (vector<Read*>::iterator r = reads.begin(); r != reads.end(); ++r) {
	Read& read = **r;
	read.buildMutualKmerHits();
	if (kmerCount > 0) {
	    read.alignWithCandidates(kmerCount, useSSW);
	} else {
	    read.alignWithReads(reads, useSSW);
	}
    }

    // now determine a traversal
    // for now, just describe ML set of paths following SW

    if (dotOutput) {
	cout << "digraph G {" << endl;
	for (vector<Read*>::iterator r = reads.begin(); r != reads.end(); ++r) {
	    (*r)->dotFormat(cout, dotPerBaseSWMin, dotSWMin);
	}
	cout << "}" << endl;
    } else {
	for (vector<Read*>::iterator r = reads.begin(); r != reads.end(); ++r) {
	    cout << **r << endl;
	}
    }

    reader.Close();
    if (!suppress_output)
        writer.Close();

    return 0;

}

