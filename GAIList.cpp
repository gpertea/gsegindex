//=============================================================================
//Read .BED datasets, and then find all overlaps
//by Jianglin Feng  09/05/2018
//Decomposition & simplication: 11/26/2018
//Radix sorting and one-pass loading based on lh3's cgranges: 6/20/2019
//-----------------------------------------------------------------------------
#include "GAIList.h"
#include "GRadixSorter.hh"
//#include "KRadixSorter.hh"

/*
// Radix Sort
#define RS_MIN_SIZE 64
#define RS_MAX_BITS 8
#define RegDataKEY(r) ((r).start)

struct RdxSortBucket {
	AIRegData *b, *e;
};

void rdx_insertsort(AIRegData *beg, AIRegData *end) {
	AIRegData *i;
	for (i = beg + 1; i < end; ++i)
		if (RegDataKEY(*i) < RegDataKEY(*(i - 1))) {
			AIRegData *j, tmp = *i;
			for (j = i; j > beg && RegDataKEY(tmp) < RegDataKEY(*(j-1)); --j)
				*j = *(j - 1);
			*j = tmp;
		}
}
void rdx_sort(AIRegData *beg, AIRegData *end, int n_bits, int s) {
	AIRegData *i;
	int size = 1<<n_bits, m = size - 1;
	RdxSortBucket *k, b[1<<RS_MAX_BITS], *be = b + size;
	assert(n_bits <= RS_MAX_BITS);
	for (k = b; k != be; ++k) k->b = k->e = beg;
	for (i = beg; i != end; ++i) ++b[RegDataKEY(*i)>>s&m].e;
	for (k = b + 1; k != be; ++k)
		k->e += (k-1)->e - beg, k->b = (k-1)->e;
	for (k = b; k != be;) {
		if (k->b != k->e) {
			RdxSortBucket *l;
			if ((l = b + (RegDataKEY(*k->b)>>s&m)) != k) {
				AIRegData tmp = *k->b, swap;
				do {
					swap = tmp; tmp = *l->b; *l->b++ = swap;
					l = b + (RegDataKEY(tmp)>>s&m);
				} while (l != k);
				*k->b++ = tmp;
			} else ++k->b;
		} else ++k;
	}
	for (b->b = beg, k = b + 1; k != be; ++k) k->b = (k-1)->e;
	if (s) {
		s = s > n_bits? s - n_bits : 0;
		for (k = b; k != be; ++k)
			if (k->e - k->b > RS_MIN_SIZE) rdx_sort(k->b, k->e, n_bits, s);
			else if (k->e - k->b > 1) rdx_insertsort(k->b, k->e);
	}
}
void radix_sort_regs(AIRegData *beg, AIRegData *end) {
	if (end - beg <= RS_MIN_SIZE) rdx_insertsort(beg, end);
	else rdx_sort(beg, end, RS_MAX_BITS, (sizeof(int32_t) - 1) * RS_MAX_BITS);
}
*/

int32_t bSearch(AIRegData* As, int32_t idxS, int32_t idxE, uint32_t qe) {   //find tE: index of the first item satisfying .s<qe from right
    int tE=-1, tL=idxS, tR=idxE-1, tM, d;
    uint32_t v;
    AIRegData* p=As;
    while((d=tR-tL)>1) {
    	tM = tL + (d>>1);
    	v=p[tM].start;
        if(v >= qe)
            tR = tM-1;
        else
            tL = tM;
    }
    if(p[tR].start < qe)
        tE = tR;
    else if(p[tL].start < qe)
        tE = tL;
    return tE;
}

void GAIList::init() {
	ctghash=new GHashMap<const char*, int32_t>();
	ctghash->resize(64);
	nctg = 0;
	mctg = 32;
	GMALLOC(ctg, mctg*sizeof(AICtgData));
	h_count=0;
	h_cap=1000000;
	GMALLOC(hits, h_cap*sizeof(uint32_t));
}


void GAIList::destroy() {
	for (int i = 0; i < nctg; ++i){
		free(ctg[i].name);
		free(ctg[i].glist);
		free(ctg[i].maxE);
	}
	free(ctg);
	if (hits) free(hits);
	delete ctghash;
}

#define CALLOC(type, len) ((type*)calloc((len), sizeof(type)))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))
#define EXPAND(a, m) { (m) = (m)? (m) + ((m)>>1) : 16; REALLOC((a), (m)); }

void GAIList::add(const char *chr, uint32_t s, uint32_t e, uint32_t payload) {
	if(s > e) return;
	bool cnew=false;
	uint64_t hidx=ctghash->addIfNew(chr, nctg, cnew);
	AICtgData *q;
	if (cnew) { //new contig
		if (nctg == mctg)
			{ EXPAND(ctg, mctg); }
		q = &ctg[nctg++];
		q->name=strdup(chr);
		q->nr=0; q->mr=64;
		GMALLOC( q->glist, (q->mr * sizeof(AIRegData)) );
		ctghash->setKey(hidx, q->name);
	} else {
		q = &ctg[ctghash->getValue(hidx)];
	}

	if (q->nr == q->mr)
		{ EXPAND(q->glist, q->mr); }
	AIRegData *p = &q->glist[q->nr++];
	p->start = s;
	p->end   = e;
	p->didx = payload;
}

//-------------------------------------------------------------------------------

void GAIList::loadBED(const char* fn) {
	gzFile fp;
	uint32_t k = 0;
	if ((fp=gzopen(fn, "r"))) {
		GFStream<gzFile, int (*)(gzFile, voidp, unsigned int)> fs(fp, gzread);
		Gcstr line;
		while (fs.getUntil(fs.SEP_LINE, line)>=0) {
			if (line.len()==0) continue;
			char *ctg;
			uint32_t st, en;
			ctg = parse_bed(line(), &st, &en);
			if (ctg) this->add(ctg, st, en, k++);

		}
	} else GError("Error: failed to open file %s\n", fn);
	gzclose(fp);
}

void GAIList::build(int cLen) {
	int cLen1=cLen/2, nr, minL = GMAX(64, cLen);
	cLen += cLen1;
	int lenT, len, iter, i, j, k, k0, t;
	GRadixSorter<AIRegData, uint32_t, & AIRegData::start> rdx;
	//KRadixSorter<AIRegData, uint32_t, & AIRegData::start> rdx;
	for(i=0; i<nctg; i++){
		//1. Decomposition
		AICtgData *p    = &ctg[i];
		AIRegData *L1 = p->glist;  //L1: to be rebuilt
		nr = p->nr;
		//radix_sort_intv(L1, L1+nr);
		//radix_sort_regs(L1, L1+nr);
		//radix_sort(L1, nr);
		rdx.sort(L1, nr);
		if(nr<=minL){
			p->nc = 1, p->lenC[0] = nr, p->idxC[0] = 0;
		}
		else{
			AIRegData *L0 = (AIRegData *)malloc(nr*sizeof(AIRegData)); 	//L0: serve as input list
			AIRegData *L2 = (AIRegData *)malloc(nr*sizeof(AIRegData));   //L2: extracted list
			//----------------------------------------
			AIRegData *D0 = (AIRegData *)malloc(nr*sizeof(AIRegData)); 	//D0:
			int32_t *di = (int32_t*)malloc(nr*sizeof(int32_t));	//int64_t?
			//----------------------------------------
			memcpy(L0, L1, nr*sizeof(AIRegData));
			iter = 0;	k = 0;	k0 = 0;
			lenT = nr;
			while(iter<MAXC && lenT>minL){
				//setup di---------------------------
				for(j=0;j<lenT;j++){				//L0:{.start= end, .end=idx, .value=idx1}
					D0[j].start = L0[j].end;
					D0[j].end = j;
				}
				//radix_sort_intv(D0, D0+lenT);
				//radix_sort_regs(D0, D0+lenT);
				//radix_sort6(D0,lenT);
				rdx.sort(D0, lenT);
				for(j=0;j<lenT;j++){				//assign i=29 to L0[i].end=2
					t = D0[j].end;
					di[t] = j-t;					//>0 indicate containment
				}
				//-----------------------------------
				len = 0;
				for(t=0;t<lenT-cLen;t++){
					if(di[t]>cLen)
						memcpy(&L2[len++], &L0[t], sizeof(AIRegData));
					else
						memcpy(&L1[k++], &L0[t], sizeof(AIRegData));
				}
				memcpy(&L1[k], &L0[lenT-cLen], cLen*sizeof(AIRegData));
				k += cLen, lenT = len;
				p->idxC[iter] = k0;
				p->lenC[iter] = k-k0;
				k0 = k, iter++;
				if(lenT<=minL || iter==MAXC-2){			//exit: add L2 to the end
					if(lenT>0){
						memcpy(&L1[k], L2, lenT*sizeof(AIRegData));
						p->idxC[iter] = k;
						p->lenC[iter] = lenT;
						iter++;
						lenT = 0;						//exit!
					}
					p->nc = iter;
				}
				else memcpy(L0, L2, lenT*sizeof(AIRegData));
			}
			free(L2),free(L0), free(D0), free(di);
		}
		//2. Augmentation
		p->maxE = (uint32_t*)malloc(nr*sizeof(uint32_t));
		for(j=0; j<p->nc; j++){
			k0 = p->idxC[j];
			k = k0 + p->lenC[j];
			uint32_t tt = L1[k0].end;
			p->maxE[k0]=tt;
			for(t=k0+1; t<k; t++){
				if(L1[t].end > tt) tt = L1[t].end;
				p->maxE[t] = tt;
			}
		}
	}
}

inline int32_t GAIList::getCtg(const char *chr) {
	int32_t* fi=this->ctghash->Find(chr);
	return fi ? *fi : -1;
}

uint32_t GAIList::query(char *chr, uint32_t qs, uint32_t qe) {
    //interestingly enough having a local copy of the hits-related fields
    // leads to better speed optimization! (likely due to using registers?)
    uint32_t nr = 0, m = h_cap, newc;
    uint32_t* r = hits;
    int32_t gid = getCtg(chr);
    if(gid>=nctg || gid<0)return 0;
    AICtgData *p = &ctg[gid];
    for(int k=0; k<p->nc; k++) { //search each component
        int32_t cs = p->idxC[k];
        int32_t ce = cs + p->lenC[k];
        int32_t t;
        if(p->lenC[k]>15){ //use binary search
            t = bSearch(p->glist, cs, ce, qe); 	//rs<qe: inline not better
            if(t>=cs){
            	newc=nr+t-cs;
		        if(newc>=m){
		        	m = newc + 1024;
		        	r = (uint32_t*)realloc(r, m*sizeof(uint32_t));
		        }
		        while(t>=cs && p->maxE[t]>qs){
		            if(p->glist[t].end>qs)
		                r[nr++] = t;
		            t--;
		        }
            }
        }
        else { //use linear search
        	newc=nr+ce-cs;
        	if(newc>=m){
        		m = newc + 1024;
        		r = (uint32_t*)realloc(r, m*sizeof(uint32_t));
        	}
            for(t=cs; t<ce; t++)
                if(p->glist[t].start<qe && p->glist[t].end>qs)
                    r[nr++] = t;
        }
    }
    //update back the hits-related fields:
    hits = r, h_cap=m, h_count=nr;
    return nr;
 }

