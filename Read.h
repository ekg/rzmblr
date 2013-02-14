#include <vector>
#include <map>
#include <string>
#include "ssw_cpp.h"
#include "smithwaterman/SmithWatermanGotoh.h"
#include "api/BamAlignment.h"
#include <sstream>
#include <algorithm>


using namespace std;
using namespace BamTools;

// utility functions
ostream& operator<<(ostream& out, vector<pair<int, char> >& cigar);
vector<pair<int, char> > splitCigar(const string& cigarStr);

class Read {

    friend ostream& operator<<(ostream& o, Read& r);

    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;

public:
    string sequence;
    string name;
    BamAlignment* bamAlignment;
    int kmerLength;
    map<string, vector<Read*> >* kmersToReads;
    map<Read*, int> mutualKmerHits;
    map<int, vector<Read*> > mutualKmerHitCounts; // sorted mutualKmerHits by count
    vector<map<string, vector<Read*> >::iterator> kmers;
    map<Read*, StripedSmithWaterman::Alignment> alignments;
    map<int, vector<Read*> > swScores;
    map<float, vector<Read*> > swScoresPerBase;
    void update();
    bool toRemove;
    void kmersOf(void);
    void buildMutualKmerHits(void);
    void alignWithCandidates(int minKmerMatchCount, bool useSSW = true);
    void alignWithReads(vector<Read*>& reads, bool useSSW = true);
    void dotFormat(ostream& out, float minPerBaseSW = 1.0, int swScoreCut = 0);

    Read(string& s,
	 map<string, vector<Read*> >* m,
	 int kl = 11,
	 BamAlignment* al = NULL,
	 int ms = 2,
	 int mm = 2,
	 int go = 3,
	 int ge = 1)
	: sequence(s)
	, bamAlignment(al)
	, kmerLength(kl)
	, kmersToReads(m)
	, toRemove(false) // for cleanup
    {
	if (bamAlignment != NULL) {
	    name = bamAlignment->Name;
	    if (bamAlignment->IsFirstMate()) {
		name += ".1";
	    } else {
		name += ".2";
	    }
	    if (bamAlignment->IsMapped()) {
		name += ".m";
	    } else {
		name += ".u";
	    }
	} else {
	    name = "unknown";
	}
	aligner.SetMismatchPenalty(ms, mm);
	aligner.SetGapPenalty(go, ge);
	kmersOf();
    }


};
