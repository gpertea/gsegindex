//=====================================================================================
//Common structs, parameters, functions
//by Jianglin Feng  09/5/2018
//-------------------------------------------------------------------------------------

#ifndef __GAILIST_H__
#define __GAILIST_H__
//-------------------------------------------------------------------------------------
#include "GHashMap.hh"

// max number of components:

//-------------------------------------------------------------------------------------

template <typename REC> class GAIList {
  protected:
	static const int MAXC=10;
	struct AIRegData {
	    uint32_t start;   //region start: 0-based
	    uint32_t end;     //region end: not inclusive
	    REC data;      //  this could be an index-1 into a storage array for data
	                   //  associated with each interval (e.g. GFF records)
	};
	struct AICtgData{
		char *name;  //name of contig
		int64_t nr, mr;  //number/capacity of glist
		AIRegData *glist;   //regions data
		int nc, lenC[MAXC], idxC[MAXC];  //components
		uint32_t *maxE;                  //augmentation
	};
	//AICtgData *ctglst; // list of contigs (of size nctg)
	AICtgData *ctg;            // list of contigs (of size nctg)
	int32_t nctg, mctg;   // count and max number of contigs
	//void *hc;              // dict for converting contig names to int
	GHashMap<const char*, int32_t>* ctghash;
	// hits buffer for the last query()
	uint32_t* hits;
	uint32_t h_cap;
	uint32_t h_count;
    //--methods
	static int32_t bSearch(AIRegData* As, int32_t idxS, int32_t idxE, uint32_t qe);
	int32_t getCtg(const char *chr) {
		int32_t* fi=this->ctghash->Find(chr);
		return fi ? *fi : -1;
	}
	void init() {
		ctghash=new GHashMap<const char*, int32_t>();
		ctghash->resize(64);
		nctg = 0;
		mctg = 32;
		GMALLOC(ctg, mctg*sizeof(AICtgData));
		h_count=0;
		h_cap=1000000;
		GMALLOC(hits, h_cap*sizeof(uint32_t));
	}
	void destroy() {
		for (int i = 0; i < nctg; ++i){
			free(ctg[i].name);
			free(ctg[i].glist);
			free(ctg[i].maxE);
		}
		free(ctg);
		if (hits) free(hits);
		delete ctghash;
	}
  public:
	void add(const char *chr, uint32_t s, uint32_t e, REC payload);
	void build(int cLen); //ailist_construct
	uint32_t query(const char *chr, uint32_t qs, uint32_t qe);
	GAIList() { this->init(); }
	~GAIList() { this->destroy(); }
};

#endif
