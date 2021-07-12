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
//#include "GHashMap.hh"
KSTREAM_INIT(gzFile, gzread, 0x10000)

//#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
//#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAXC 10							//max number of components


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

struct AIList {
	ctg_t *ctg;        					// list of contigs (of size _n_ctg_)
	int32_t nctg, mctg; 				// number and max number of contigs
	//void *hc;             				// dict for converting contig names to int
	GHashMap<const char*, int32_t>* ctghash;
};
/*
struct AICtgData{
	char *ctg=NULL;            //name of the contig
	uint32_t nr=0, maxr=64;    //number of regions
	AIData *reglist=NULL;        //regions data
	//GDynArray<AIData> reglist;
	uint32_t nc=0;               //number of components
	uint lenC[MAXC], idxC[MAXC]; //components
	uint32_t *maxE=NULL;      //augmentation
	AICtgData(const char* chr=NULL) {
		init(chr);
	}
	void init(const char* chr) {
		nr=0;maxr=64;
		nc=0;
		ctg=strdup(chr);
		GMALLOC(reglist, maxr * sizeof(AIData));
	}
	AIData& add() {
		if (nr==maxr) {
			maxr = maxr ? maxr + (maxr>>1) : 16;
			GREALLOC(reglist, maxr * sizeof(AIData));
		}
		nr++;
		return reglist[nr-1];
	}

	void destroy() {
		free(ctg);
		free(maxE);
		free(reglist);
	}
};

*/

struct GAIList{
	//AICtgData *ctglst; // list of contigs (of size nctg)
	ctg_t *ctg;            // list of contigs (of size nctg)
	uint32_t nctg, mctg;   // count and max number of contigs
	//void *hc;              // dict for converting contig names to int
	GHashMap<const char*, int32_t>* ctghash;

	void init();
	void loadBED(const char* fn);
	void add(const char *chr, uint32_t s, uint32_t e, uint32_t payload);
	void build(int cLen); //ailist_construct
	void query(char *chr, uint32_t qs, uint32_t qe, uint32_t *mr, uint32_t **ir);
	int32_t getCtg(const char *chr);
	void destroy();
	GAIList() { init(); }
	~GAIList() { destroy(); }
};
/*
struct GAIList {
	//GPVec<AICtgData> ctglst;        // list of contigs (of size _n_ctg_)
	GDynArray<AICtgData> ctgs;
	//uint32_t nctg, mctg;   // number and max number of contigs
	void *hc=NULL;              // hashmap for converting contig names to int
	GAIList(uint max=32); //init
	void loadBED(const char* fn); //readBED
	void add(const char *chr, uint32_t s, uint32_t e, uint32_t payload); //ailist_add
	void build(int cLen); //ailist_construct
	int32_t getCtg(const char* chr);
	// returns an array of indexes mr
	uint32_t query(char *chr, uint32_t qs, uint32_t qe, GDynArray<uint32_t>& hits); //ailist_query
    ~GAIList();
};
*/
//-------------------------------------------------------------------------------------

//Initialize AIList
AIList *gailist_init(void);

//read .BED file
AIList* greadBED(const char* fn);

//Add a AIData interval
void gailist_add(AIList *ail, const char *chr, uint32_t s, uint32_t e, uint32_t v);

//Construct ailist: decomposition and augmentation
void gailist_construct(AIList *ail, int cLen);
//void ailist_construct0(AIList *ail, int cLen);

//Get chr index
//int32_t get_ctg(const void* hp, const char *chr);
//int32_t get_ctg(const AIList *ail, const char *chr);
int32_t get_ctg(const AIList *ail, const char *chr);

//Binary search
uint32_t bSearch(AIData* As, uint32_t idxS, uint32_t idxE, uint32_t qe);

//Query ailist intervals
uint32_t gailist_query(AIList *ail, char *chr, uint32_t qs, uint32_t qe, uint32_t *mr, uint32_t **ir);

//Free ailist data
void gailist_destroy(AIList *ail);
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

/*********************
 * Convenient macros *
 *********************/

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define CALLOC(type, len) ((type*)calloc((len), sizeof(type)))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

#define EXPAND(a, m) do { \
		(m) = (m)? (m) + ((m)>>1) : 16; \
		REALLOC((a), (m)); \
	}while (0)

#endif
