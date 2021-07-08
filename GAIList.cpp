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

GAIList *ailist_init(void)
{
	GAIList *ail = NULL;
	GMALLOC(ail,  1 * sizeof(GAIList));
	ail->hc = kh_init(str);
	ail->nctg = 0;
	ail->mctg = 32; //preallocate mem for 32 contigs
	GMALLOC(ail->ctglst, ail->mctg *sizeof(AICtgData));
	return ail;
}

void ailist_destroy(GAIList *ail)
{
	uint32_t i;
	if (ail == 0) return;
	for (i = 0; i < ail->nctg; ++i){
		free(ail->ctglst[i].ctg);
		free(ail->ctglst[i].glist);
		free(ail->ctglst[i].maxE);
	}
	free(ail->ctglst);
	kh_destroy(str, (strhash_t*)ail->hc);
	free(ail);
}

void ailist_add(GAIList *ail, const char *chr, uint32_t s, uint32_t e) //, int32_t v)
{
	if(s > e)return;
	int absent;
	khint_t k;
	strhash_t *h = (strhash_t*)ail->hc;
	k = kh_put(str, h, chr, &absent);
	if (absent) {
		if (ail->nctg == ail->mctg)
			EXPAND(ail->ctglst, ail->mctg);
		kh_val(h, k) = ail->nctg;
		AICtgData *p = &ail->ctglst[ail->nctg++];
		p->ctg = strdup(chr);
		p->nr=0;	p->maxr=64;
		GMALLOC(p->glist, p->maxr * sizeof(AIData));
		kh_key(h, k) = p->ctg;
	}
	int32_t kk = kh_val(h, k);
	AICtgData *q = &ail->ctglst[kk];
	if (q->nr == q->maxr)
		EXPAND(q->glist, q->maxr);
	AIData *p = &q->glist[q->nr++];
	p->start = s;
	p->end   = e;
	return;
}

//-------------------------------------------------------------------------------

void GAIList::loadBED(const char* fn) {
	gzFile fp;
	int32_t k = 0;
	if ((fp=gzopen(fn, "r"))) {
		GFStream<gzFile, int (*)(gzFile, voidp, unsigned int)> fs(fp, gzread);
		Gcstr line;
		while (fs.getUntil(fs.SEP_LINE, line)>=0) {
			if (line.len()==0) continue;
			char *ctg;
			int32_t st, en;
			ctg = parse_bed(line(), &st, &en);
			if (ctg) this->add(ctg, st, en, k++);

		}
	} else GError("Error: failed to open file %s\n", fn);
	gzclose(fp);
}

void GAIList::add(const char *chr, uint32_t s, uint32_t e, uint32_t payload) {
	if(s > e)return;
	int absent;
	khint_t k;
	strhash_t *h = (strhash_t*)(this->hc);
	k = kh_put(str, h, chr, &absent);
	if (absent) { //add new AICtgData
		//kh_val(h, k) = ctglst.Count();
		AICtgData* newctg=new AICtgData(chr);
		uint ctgidx=ctglst.Add(newctg);
		kh_val(h, k) = ctgidx;
		kh_key(h, k) = newctg->ctg;
	}
	int32_t kk = kh_val(h, k);
	AICtgData* q = ctglst[kk];
	/*
	if (q.nr == q.maxr)
		EXPAND(q.glist, q.maxr);
	AIData *p = &q->glist[q->nr++];
	p->start = s;
	p->end   = e;
	*/
	q->reglist.Add({s, e, payload});

}


// TFunc for GFStream:
//typedef int (*ft_gzread)(gzFile, voidp, unsigned int);
GAIList* readBED(const char* fn)
{   //faster than strtok()
	gzFile fp;
	GAIList *ail;
	//int32_t k = 0;
	ail = ailist_init();
	if ((fp=gzopen(fn, "r"))) {
		GFStream<gzFile, int (*)(gzFile, voidp, unsigned int)> fs(fp, gzread);
		Gcstr line;
		while (fs.getUntil(fs.SEP_LINE, line)>=0) {
			if (line.len()==0) continue;
			char *ctg;
			int32_t st, en;
			ctg = parse_bed(line(), &st, &en);
			if (ctg) ailist_add(ail, ctg, st, en);//, k++);

		}
	} else GError("Error: failed to open file %s\n", fn);
	gzclose(fp);
	return ail;
}


void ailist_construct(GAIList *ail, int cLen) {
	int cLen1=cLen/2;
	uint32_t minL = MAX(64, cLen);
	uint32_t nr;
	cLen += cLen1;
	uint32_t lenT, len, iter, j, i, t, k, k0;
	for(i=0; i<ail->nctg; i++){ //for each chromosome
		//1. Decomposition
		AICtgData *ctgdta    = &ail->ctglst[i];
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

void GAIList::build(int cLen) {
	int cLen1=cLen/2;
	uint32_t minL = MAX(64, cLen);
	uint32_t nr;
	cLen += cLen1;
	uint32_t lenT, len, iter, j, t, k, k0;
	for(uint i=0; i<ctglst.Count(); i++){ //for each chromosome
		//1. Decomposition
		AICtgData& ctgdta = ctglst[i];
		AIData* L1 = ctgdta.reglist(); //L1: to be rebuilt
		nr = ctgdta.reglist.Count();
		radix_sort_intv(L1, L1+nr);
		if (nr<=minL) {
			ctgdta.nc = 1;
			ctgdta.lenC[0] = nr;
			ctgdta.idxC[0] = 0;
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
				ctgdta.idxC[iter] = k0;
				ctgdta.lenC[iter] = k-k0;
				k0 = k, iter++;
				if (lenT<=minL || iter==MAXC-2){			//exit: add L2 to the end
					if(lenT>0){
						memcpy(&L1[k], L2, lenT*sizeof(AIData));
						ctgdta.idxC[iter] = k;
						ctgdta.lenC[iter] = lenT;
						iter++;
						lenT = 0;						//exit!
					}
					ctgdta.nc = iter;
				}
				else memcpy(L0, L2, lenT*sizeof(AIData));
			}
			free(L2),free(L0), free(D0), free(di);
		}
		//2. Augmentation
		GMALLOC(ctgdta.maxE,  nr*sizeof(uint32_t));
		for(j=0; j<ctgdta.nc; j++){
			k0 = ctgdta.idxC[j];
			k = k0 + ctgdta.lenC[j];
			uint32_t tt = L1[k0].end;
			ctgdta.maxE[k0]=tt;
			for(t=k0+1; t<k; t++){
				if(L1[t].end > tt) tt = L1[t].end;
				ctgdta.maxE[t] = tt;
			}
		}
	}
}


int32_t GAIList::get_ctg(const char *chr) {
	khint_t k;
	strhash_t *h = (strhash_t*)hc;
	k = kh_get(str, h, chr);
	return k == kh_end(h) ? -1 : kh_val(h, k);
}


int32_t get_ctg(const GAIList *ail, const char *chr)
{
	khint_t k;
	strhash_t *h = (strhash_t*)ail->hc;
	k = kh_get(str, h, chr);
	return k == kh_end(h)? -1 : kh_val(h, k);
}


uint32_t GAIList::query(char *chr, uint32_t qs, uint32_t qe, uint32_t *mr, uint32_t **ir) {
    uint32_t nr = 0, m = *mr, *r = *ir;
    int32_t gid = this->get_ctg(chr);

    if (gid<0 || (uint32_t)gid >= ctglst.Count() )
    	return 0; //no such contig
    AICtgData &p = ctglst[gid];
    for(uint k=0; k<p->nc; k++){     //search each component
        int32_t cs = p->idxC[k];
        int32_t ce = cs + p->lenC[k];
        int32_t t;
        if(p->lenC[k]>15){
            t = bSearch(p->glist, cs, ce, qe);  //rs<qe: inline not better
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

uint32_t ailist_query(GAIList *ail, char *chr, uint32_t qs, uint32_t qe, uint32_t *mr, uint32_t **ir)
{
    uint32_t nr = 0, m = *mr, *r = *ir;
    int32_t gid = get_ctg(ail, chr);

    if (gid<0 || (uint32_t)gid >= ail->nctg )
    	return 0; //no such contig
    AICtgData *p = &ail->ctglst[gid];
    for(uint k=0; k<p->nc; k++){     //search each component
        int32_t cs = p->idxC[k];
        int32_t ce = cs + p->lenC[k];
        int32_t t;
        if(p->lenC[k]>15){
            t = bSearch(p->glist, cs, ce, qe);  //rs<qe: inline not better
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

