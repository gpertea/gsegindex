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
struct AIData {
    uint32_t start;   //region start: 0-based
    uint32_t end;     //region end: not inclusive
    uint32_t didx; // this could be an index-1 into a storage array for data
                  //      associated with each interval (e.g. GFF records)
};

typedef struct{
	char *name;    						//name of the contig
	int64_t nr, mr;						//number of regions
	AIData *glist;						//regions data
	int nc, lenC[MAXC], idxC[MAXC];		//components
	uint32_t *maxE;						//augmentation
} ctg_t;

struct GAIList{
	//AICtgData *ctglst; // list of contigs (of size nctg)
	ctg_t *ctg;            // list of contigs (of size nctg)
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

//-------------------------------------------------------------------------------------
//The following section taken from Dr Heng Li's cgranges
// (https://github.com/lh3/cgranges)

//KSTREAM_INIT(gzFile, gzread, 0x10000)
/**************
 * Radix sort *
 **************/
#define RS_MIN_SIZE 64
#define RS_MAX_BITS 8

#define KRADIX_SORT_INIT(name, rstype_t, rskey, sizeof_key) \
	typedef struct { \
		rstype_t *b, *e; \
	} rsbucket_##name##_t; \
	void rs_insertsort_##name(rstype_t *beg, rstype_t *end) \
	{ \
		rstype_t *i; \
		for (i = beg + 1; i < end; ++i) \
			if (rskey(*i) < rskey(*(i - 1))) { \
				rstype_t *j, tmp = *i; \
				for (j = i; j > beg && rskey(tmp) < rskey(*(j-1)); --j) \
					*j = *(j - 1); \
				*j = tmp; \
			} \
	} \
	void rs_sort_##name(rstype_t *beg, rstype_t *end, int n_bits, int s) \
	{ \
		rstype_t *i; \
		int size = 1<<n_bits, m = size - 1; \
		rsbucket_##name##_t *k, b[1<<RS_MAX_BITS], *be = b + size; \
		assert(n_bits <= RS_MAX_BITS); \
		for (k = b; k != be; ++k) k->b = k->e = beg; \
		for (i = beg; i != end; ++i) ++b[rskey(*i)>>s&m].e; \
		for (k = b + 1; k != be; ++k) \
			k->e += (k-1)->e - beg, k->b = (k-1)->e; \
		for (k = b; k != be;) { \
			if (k->b != k->e) { \
				rsbucket_##name##_t *l; \
				if ((l = b + (rskey(*k->b)>>s&m)) != k) { \
					rstype_t tmp = *k->b, swap; \
					do { \
						swap = tmp; tmp = *l->b; *l->b++ = swap; \
						l = b + (rskey(tmp)>>s&m); \
					} while (l != k); \
					*k->b++ = tmp; \
				} else ++k->b; \
			} else ++k; \
		} \
		for (b->b = beg, k = b + 1; k != be; ++k) k->b = (k-1)->e; \
		if (s) { \
			s = s > n_bits? s - n_bits : 0; \
			for (k = b; k != be; ++k) \
				if (k->e - k->b > RS_MIN_SIZE) rs_sort_##name(k->b, k->e, n_bits, s); \
				else if (k->e - k->b > 1) rs_insertsort_##name(k->b, k->e); \
		} \
	} \
	void radix_sort_##name(rstype_t *beg, rstype_t *end) \
	{ \
		if (end - beg <= RS_MIN_SIZE) rs_insertsort_##name(beg, end); \
		else rs_sort_##name(beg, end, RS_MAX_BITS, (sizeof_key - 1) * RS_MAX_BITS); \
	}


#endif
