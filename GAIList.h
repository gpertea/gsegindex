//=====================================================================================
//Common structs, parameters, functions
//by Jianglin Feng  09/5/2018
//-------------------------------------------------------------------------------------

#ifndef __GAILIST_H__
#define __GAILIST_H__
//-------------------------------------------------------------------------------------
#include "iutil.h"
#include "GVec.hh"
#include "GHashMap.hh"

// max number of components:
#define MAXC 10

//-------------------------------------------------------------------------------------
struct AIRegData {
    uint32_t start;   //region start: 0-based
    uint32_t end;     //region end: not inclusive
    uint32_t didx; // this could be an index-1 into a storage array for data
                  //      associated with each interval (e.g. GFF records)
};

struct AICtgData{
	char *name;  //name of contig
	int64_t nr, mr;  //number/capacity of glist
	AIRegData *glist;   //regions data
	int nc, lenC[MAXC], idxC[MAXC];  //components
	uint32_t *maxE;                  //augmentation
};

struct GAIList{
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
	void init();
	void loadBED(const char* fn);
	void add(const char *chr, uint32_t s, uint32_t e, uint32_t payload);
	void build(int cLen); //ailist_construct
	uint32_t query(char *chr, uint32_t qs, uint32_t qe);
	int32_t getCtg(const char *chr);
	void destroy();
	GAIList() { init(); }
	~GAIList() { destroy(); }
};

#endif
