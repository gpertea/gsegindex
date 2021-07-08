//=============================================================================
//Read .BED datasets, and then find all overlaps
//by Jianglin Feng  09/05/2018
//Decomposition & simplication: 11/26/2018
//Radix sorting and one-pass loading based on lh3's cgranges: 6/20/2019
//-----------------------------------------------------------------------------
#include "GAIList.h"

#define AIData_key(r) ((r).start)
KRADIX_SORT_INIT(intv, AIData, AIData_key, 4)

KHASH_MAP_INIT_STR(str, int32_t)
typedef khash_t(str) strhash_t;

uint32_t bSearch(AIData* As, uint32_t idxS, uint32_t idxE, uint32_t qe) {   //find tE: index of the first item satisfying .s<qe from right
    int tL=idxS, tR=idxE-1, tM, tE=-1;
    if(As[tR].start < qe)
        return tR;
    else if(As[tL].start >= qe)
        return -1;
    while(tL<tR-1){
        tM = (tL+tR)/2;
        if(As[tM].start >= qe)
            tR = tM-1;
        else
            tL = tM;
    }
    if(As[tR].start < qe)
        tE = tR;
    else if(As[tL].start < qe)
        tE = tL;
    return tE;
}

AIList *ailist_init(void)
{
	AIList *ail = NULL;
	GMALLOC(ail,  1 * sizeof(AIList));
	ail->hc = kh_init(str);
	ail->nctg = 0;
	ail->mctg = 32;
	GMALLOC(ail->ctg, ail->mctg *sizeof(AICtgData));
	return ail;
}

void ailist_destroy(AIList *ail)
{
	uint32_t i;
	if (ail == 0) return;
	for (i = 0; i < ail->nctg; ++i){
		free(ail->ctg[i].ctg);
		free(ail->ctg[i].glist);
		free(ail->ctg[i].maxE);
	}
	free(ail->ctg);
	kh_destroy(str, (strhash_t*)ail->hc);
	free(ail);
}

void ailist_add(AIList *ail, const char *chr, uint32_t s, uint32_t e, int32_t v)
{
	if(s > e)return;
	int absent;
	khint_t k;
	strhash_t *h = (strhash_t*)ail->hc;
	k = kh_put(str, h, chr, &absent);
	if (absent) {
		if (ail->nctg == ail->mctg)
			EXPAND(ail->ctg, ail->mctg);
		kh_val(h, k) = ail->nctg;
		AICtgData *p = &ail->ctg[ail->nctg++];
		p->ctg = strdup(chr);
		p->nr=0;	p->maxr=64;
		GMALLOC(p->glist, p->maxr * sizeof(AIData));
		kh_key(h, k) = p->ctg;
	}
	int32_t kk = kh_val(h, k);
	AICtgData *q = &ail->ctg[kk];
	if (q->nr == q->maxr)
		EXPAND(q->glist, q->maxr);
	AIData *p = &q->glist[q->nr++];
	p->start = s;
	p->end   = e;
	return;
}

//-------------------------------------------------------------------------------

// TFunc for GFStream:
//typedef int (*ft_gzread)(gzFile, voidp, unsigned int);
AIList* readBED(const char* fn)
{   //faster than strtok()
	gzFile fp;
	AIList *ail;
	/*
	kstream_t *ks;
	kstring_t str = {0,0,0};
	*/
	int32_t k = 0;

	ail = ailist_init();
	/*
	if ((fp = gzopen(fn, "r")) == 0)
		return 0;
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		char *ctg;
		int32_t st, en;
		ctg = parse_bed(str.s, &st, &en);
		if (ctg) ailist_add(ail, ctg, st, en, k++);
	}
	free(str.s);
	ks_destroy(ks);
	*/
	if ((fp=gzopen(fn, "r"))) {
		GFStream<gzFile, int (*)(gzFile, voidp, unsigned int)> fs(fp, gzread);
		Gcstr line;
		while (fs.getUntil(fs.SEP_LINE, line)>=0) {
			if (line.len()==0) continue;
			char *ctg;
			int32_t st, en;
			ctg = parse_bed(line(), &st, &en);
			if (ctg) ailist_add(ail, ctg, st, en, k++);

		}
	}
	gzclose(fp);
	return ail;
}

void ailist_construct(AIList *ail, int cLen) {
	int cLen1=cLen/2;
	uint32_t minL = MAX(64, cLen);
	uint32_t nr;
	cLen += cLen1;
	uint32_t lenT, len, iter, j, i, t, k, k0;
	for(i=0; i<ail->nctg; i++){ //for each chromosome
		//1. Decomposition
		AICtgData *ctgdta    = &ail->ctg[i];
		AIData *L1 = ctgdta->glist;							//L1: to be rebuilt
		nr = ctgdta->nr;
		radix_sort_intv(L1, L1+nr);
		if (nr<=minL) {
			ctgdta->nc = 1;
			ctgdta->lenC[0] = nr;
			ctgdta->idxC[0] = 0;
		}
		else{
			AIData *L0 = malloc(nr*sizeof(AIData)); 	//L0: serve as input list
			AIData *L2 = malloc(nr*sizeof(AIData));   //L2: extracted list
			//----------------------------------------
			AIData *D0 = malloc(nr*sizeof(AIData)); 	//D0:
			int32_t *di = malloc(nr*sizeof(int32_t));	//int64_t?
			//----------------------------------------
			memcpy(L0, L1, nr*sizeof(AIData));
			iter = 0;	k = 0;	k0 = 0;
			lenT = nr;
			while(iter<MAXC && lenT>minL){
				//setup di---------------------------
				for(j=0;j<lenT;j++){				//L0:{.start= end, .end=idx, .value=idx1}
					D0[j].start = L0[j].end;
					D0[j].end = j;
				}
				radix_sort_intv(D0, D0+lenT);
				for(j=0;j<lenT;j++){				//assign i=29 to L0[i].end=2
					t = D0[j].end;
					di[t] = j-t;					//>0 indicate containment
				}
				//-----------------------------------
				len = 0;
				for(t=0;t<lenT-cLen;t++){
					if(di[t]>cLen)
						memcpy(&L2[len++], &L0[t], sizeof(AIData));
					else
						memcpy(&L1[k++], &L0[t], sizeof(AIData));
				}
				memcpy(&L1[k], &L0[lenT-cLen], cLen*sizeof(AIData));
				k += cLen;
				lenT = len;
				ctgdta->idxC[iter] = k0;
				ctgdta->lenC[iter] = k-k0;
				k0 = k, iter++;
				if (lenT<=minL || iter==MAXC-2){			//exit: add L2 to the end
					if(lenT>0){
						memcpy(&L1[k], L2, lenT*sizeof(AIData));
						ctgdta->idxC[iter] = k;
						ctgdta->lenC[iter] = lenT;
						iter++;
						lenT = 0;						//exit!
					}
					ctgdta->nc = iter;
				}
				else memcpy(L0, L2, lenT*sizeof(AIData));
			}
			free(L2),free(L0), free(D0), free(di);
		}
		//2. Augmentation
		GMALLOC(ctgdta->maxE,  nr*sizeof(uint32_t));
		for(j=0; j<ctgdta->nc; j++){
			k0 = ctgdta->idxC[j];
			k = k0 + ctgdta->lenC[j];
			uint32_t tt = L1[k0].end;
			ctgdta->maxE[k0]=tt;
			for(t=k0+1; t<k; t++){
				if(L1[t].end > tt) tt = L1[t].end;
				ctgdta->maxE[t] = tt;
			}
		}
	}
}

int32_t get_ctg(const AIList *ail, const char *chr)
{
	khint_t k;
	strhash_t *h = (strhash_t*)ail->hc;
	k = kh_get(str, h, chr);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

uint32_t ailist_query(AIList *ail, char *chr, uint32_t qs, uint32_t qe, uint32_t *mr, uint32_t **ir)
{
    uint32_t nr = 0, m = *mr, *r = *ir;
    int32_t gid = get_ctg(ail, chr);
    if (gid<0 || (uint32_t)gid >= ail->nctg )
    	return 0; //no such contig
    AICtgData *p = &ail->ctg[gid];
    for(uint k=0; k<p->nc; k++){					//search each component
        int32_t cs = p->idxC[k];
        int32_t ce = cs + p->lenC[k];
        int32_t t;
        if(p->lenC[k]>15){
            t = bSearch(p->glist, cs, ce, qe); 	//rs<qe: inline not better
            if(t>=cs){
		        if(nr+t-cs>=m){
		        	m = nr+t-cs + 1024;
		        	GREALLOC(r, m*sizeof(uint32_t));
		        }
		        while(t>=cs && p->maxE[t]>qs){
		            if(p->glist[t].end>qs)
		                r[nr++] = t;
		            t--;
		        }
            }
        }
        else{
        	if(nr+ce-cs>=m){
        		m = nr+ce-cs + 1024;
        		GREALLOC(r, m*sizeof(uint32_t));
        	}
            for(t=cs; t<ce; t++)
                if(p->glist[t].start<qe && p->glist[t].end>qs)
                    r[nr++] = t;
        }
    }
    *ir = r, *mr = m;
    return nr;
}

