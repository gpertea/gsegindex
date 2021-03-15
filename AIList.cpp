//=============================================================================
//Read .BED datasets, and then find all overlaps
//by Jianglin Feng  09/05/2018
//Decomposition & simplication: 11/26/2018
//Radix sorting and one-pass loading based on lh3's cgranges: 6/20/2019
//-----------------------------------------------------------------------------
#include "AIList.h"


#define gdata_t_key(r) ((r).start)
//KRADIX_SORT_INIT(intv, gdata_t, gdata_t_key, 4)
KRADIX_SORT_INIT(intv, GSeg, gdata_t_key, 4)

KHASH_MAP_INIT_STR(str, int32_t)

typedef khash_t(str) strhash_t;


//#define KRADIX_SORT_INIT(name, rstype_t, rskey, sizeof_key)
struct rsbucket_seg_t {
	int ib;
	int ie;
};

void rs_insertsort_seg(GPVec<GSeg>& lst, int ibeg, int iend) {
	//GSeg *i;
	for (int i = ibeg + 1; i < iend; ++i)
		if (lst[i]->start < lst[i-1]->start) {
			GSeg *stmp = lst[i];
			int j;
			for (j = i; j > ibeg && stmp->start < lst[j-1]->start; --j)
				lst[j] = lst[j - 1];
			lst[j] = stmp;
		}
}

void rs_sort_seg(GPVec<GSeg>& lst, int ibeg, int iend, int n_bits, int s) {
	int size = 1<<n_bits, m = size - 1;
	rsbucket_seg_t *k, b[1<<RS_MAX_BITS], *be = b + size;
	assert(n_bits <= RS_MAX_BITS);
	for (k = b; k != be; ++k) k->ib = k->ie = ibeg;
	for (int i = ibeg; i != iend; ++i) ++b[(lst[i]->start)>>s&m].ie;
	for (k = b + 1; k != be; ++k)
		k->ie += (k-1)->ie - ibeg, k->ib = (k-1)->ie;
	for (k = b; k != be;) {
		if (k->ib != k->ie) {
			rsbucket_seg_t *l;
			if ((l = b + ((lst[k->ib]->start)>>s&m)) != k) {
				int tmp = k->ib, swap;
				do {
					lst.Swap(tmp, l->ib);
					swap = tmp; tmp = l->ib; l->ib = swap;
					l->ib++;
					l = b + (lst[tmp]->start >>s&m);
				} while (l != k);
				k->ib = tmp;k->ib++;
			} else ++k->ib;
		} else ++k;
	}
	for (b->ib = ibeg, k = b + 1; k != be; ++k) k->ib = (k-1)->ie;
	if (s) {
		s = s > n_bits? s - n_bits : 0;
		for (k = b; k != be; ++k)
			if (k->ie - k->ib > RS_MIN_SIZE) rs_sort_seg(lst, k->ib, k->ie, n_bits, s);
			else if (k->ie - k->ib > 1) rs_insertsort_seg(lst, k->ib, k->ie);
	}
}

void radix_sort_seg(GPVec<GSeg>& lst, int ibeg, int iend) {
	if (iend - ibeg <= RS_MIN_SIZE) rs_insertsort_seg(lst, ibeg, iend);
	else rs_sort_seg(lst, ibeg, iend, RS_MAX_BITS, (sizeof(int) - 1) * RS_MAX_BITS);
}

char *parse_bed(char *s, int32_t *st_, int32_t *en_) {
	char *p, *q, *ctg = 0;
	int32_t i, st = -1, en = -1;
	for (i = 0, p = q = s;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
			if (i == 0) ctg = p;
			else if (i == 1) st = atol(p);
			else if (i == 2) en = atol(p);
			++i, p = q + 1;
			if (c == '\0') break;
		}
	}
	*st_ = st, *en_ = en;
	return i >= 3? ctg : 0;
}

uint32_t bSearch(GPVec<GSeg>& As, uint32_t idxS, uint32_t idxE, uint32_t qe)
{   //find tE: index of the first item satisfying .s<qe from right
    int tL=idxS, tR=idxE-1, tM, tE=-1;
    if(As[tR]->start < qe)
        return tR;
    else if(As[tL]->start >= qe)
        return -1;
    while(tL<tR-1){
        tM = tL + ((tR-tL)>>1);
        if(As[tM]->start >= qe)
            tR = tM-1;
        else
            tL = tM;
    }
    if(As[tR]->start < qe)
        tE = tR;
    else if(As[tL]->start < qe)
        tE = tL;
    return tE;
}

/* TAIList *ailist_init(void)
{
	TAIList *ail = malloc(1*sizeof(TAIList));
	ail->hc = kh_init(str);
	ail->nctg = 0;
	ail->mctg = 32;
	ail->ctg = malloc(ail->mctg*sizeof(TAILstCtg));
	return ail;
}*/

TAIList::TAIList():ctg(10, true),hc(NULL) {
	//initialization code
	hc = kh_init(str);
	//nctg = 0;
	//mctg = 32;
	//ctg = malloc(ail->mctg*sizeof(TAILstCtg));
}

TAIList::~TAIList() {
	kh_destroy(str, (strhash_t*)hc);
}

void TAIList::add(const char *chr, uint32_t s, uint32_t e, int32_t v) {
	if(s > e)return;
	int absent;
	khint_t k;
	k = kh_put(str, (strhash_t*)hc, chr, &absent);
	if (absent) { //add this contig
		/* if (ail->nctg == ail->mctg)
			EXPAND(ail->ctg, ail->mctg);*/
		//kh_val(h, k) = ail->nctg;
		TAILstCtg* nc=new TAILstCtg(chr);
		kh_val((strhash_t*)hc, k) = ctg.Add(nc);
		/* TAILstCtg *p = &ail->ctg[ail->nctg++];
		p->name = strdup(chr);
		p->nr=0;	p->mr=64;
		p->glist = malloc(p->mr*sizeof(gdata_t));
		*/
		kh_key((strhash_t*)hc, k) = nc->name;
	}
	int32_t kk = kh_val((strhash_t*)hc, k);
	TAILstCtg *q = ctg[kk];
	/* if (q->nr == q->mr) {
		EXPAND(q->glist, q->mr);
	}

	GSeg *p = q->glist[q->nr++]; //&q->glist[q->nr++];
    p->start = s;
	p->end   = e;
	*/
	GSeg *p = new GSeg(s, e);
	q->glist.Add(p);
	return;
}

//-------------------------------------------------------------------------------
void readBED(TAIList& ail, const char* fn) {
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	int32_t k = 0;
	if ((fp = gzopen(fn, "r")) == 0)
		return;
	ks = ks_init(fp);
	//ail = ailist_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		char *ctg;
		int32_t st, en;
		ctg = parse_bed(str.s, &st, &en);
		if (ctg) ail.add(ctg, st, en, k++);
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	//return ail;
}

void TAIList::build(int cLen) {
    int cLen1=cLen/2, nr, minL = GMAX(64, cLen);
    cLen += cLen1;
    int lenT, len, iter, i, j, k, k0, t;
	for(i=0; i<ctg.Count(); i++){
		//1. Decomposition
		TAILstCtg *p    = ctg[i];
		//gdata_t *L1 = p->glist;							//L1: to be rebuilt
		GPVec<GSeg>& L1 = p->glist;
		//nr 			= p->nr;
		nr=p->glist.Count();
		radix_sort_seg(p->glist, 0, nr);
		/*
		GMessage("after sorting:\n");
		for (int v=0;v<p->glist.Count();v++) {
			GMessage("[%d - %d]\n", p->glist[v]->start, p->glist[v]->end);
		}
		*/
        if(nr<=minL){
            p->nc = 1, p->lenC[0] = nr, p->idxC[0] = 0;
        }
        else{
        	/*
        	gdata_t *L0 = malloc(nr*sizeof(gdata_t)); 	//L0: serve as input list
            gdata_t *L2 = malloc(nr*sizeof(gdata_t));   //L2: extracted list
            //----------------------------------------
        	gdata_t *D0 = malloc(nr*sizeof(gdata_t)); 	//D0:
			int32_t *di = malloc(nr*sizeof(int32_t));	//int64_t?
            //----------------------------------------
            memcpy(L0, L1, nr*sizeof(gdata_t));
            */
        	GPVec<GSeg> L0(nr); L0.addNew(L1);//deep copy
        	GPVec<GSeg> L2(nr); for (int i=0;i<nr;i++) L2.Add(new GSeg());
        	GPVec<GSeg> D0(nr); //for (int i=0;i<nr;i++) D0.Add(new GSeg());
        	GVec<int> di(nr, 0);
            iter = 0;	k = 0;	k0 = 0;
            lenT = nr;
            while(iter<MAXC && lenT>minL){
            	//setup di---------------------------
		        for(j=0;j<lenT;j++){				//L0:{.start= end, .end=idx, .value=idx1}
		        	int dj=D0.Add(new GSeg(0, j));
					//D0[j].start = L0[j].end;
		        	D0[dj]->start=L0[j]->end;
					//D0[j].end = j;
				}
				radix_sort_seg(D0, 0, lenT);
				for(j=0;j<lenT;j++){				//assign i=29 to L0[i].end=2
					t = D0[j]->end;
					di[t] = j-t;					//>0 indicate containment
				}
				//-----------------------------------
                len = 0;
		        for(t=0;t<lenT-cLen;t++){
					if(di[t]>cLen) {
				        //memcpy(&L2[len++], &L0[t], sizeof(gdata_t));
						*(L2[len])=*(L0[t]);
						len++;
					}
					else {
						//memcpy(&L1[k++], &L0[t], sizeof(gdata_t));
						*(L1[k])=*(L0[t]);
						k++;
					}
				}
                //memcpy(&L1[k], &L0[lenT-cLen], cLen*sizeof(gdata_t));
		        for (int l=0;l<cLen;l++)
		        	*(L1[k+l])=*(L0[lenT-cLen+l]);
                k += cLen, lenT = len;
                p->idxC[iter] = k0;
                p->lenC[iter] = k-k0;
                k0 = k, iter++;
                if(lenT<=minL || iter==MAXC-2){			//exit: add L2 to the end
                    if(lenT>0){
                        //memcpy(&L1[k], L2, lenT*sizeof(gdata_t));
                    	for (int l=0;l<lenT;l++)
                    	   *(L1[k+l])=*(L2[l]);
                        p->idxC[iter] = k;
                        p->lenC[iter] = lenT;
                        iter++;
                        lenT = 0;						//exit!
                    }
                   	p->nc = iter;
                }
                else {
                	//memcpy(L0, L2, lenT*sizeof(gdata_t));
                	for (int l=0;l<lenT;l++)
                		*(L0[l])=*(L2[l]);
                }
            }
            //free(L2),free(L0), free(D0), free(di);
        }
        //2. Augmentation
        GMALLOC(p->maxE, nr*sizeof(uint32_t));
        for(j=0; j<p->nc; j++){
            k0 = p->idxC[j];
            k = k0 + p->lenC[j];
            uint32_t tt = L1[k0]->end;
            p->maxE[k0]=tt;
            for(t=k0+1; t<k; t++){
                if(L1[t]->end > tt) tt = L1[t]->end;
                p->maxE[t] = tt;
            }
        }
        /*
  	    GMessage("after augmentation:\n");
    	for (int v=0;v<p->glist.Count();v++) {
    		GMessage("[%d - %d]\n", p->glist[v]->start, p->glist[v]->end);
     	}
    	GMessage("%d components: \n", p->nc);
    	GMessage("maxE: (");for (int c=0;c<p->nc;c++) GMessage(" %d", p->maxE[c]);
    	GMessage(" )\n");
    	*/
	}

}
/*
void ailist_construct0(TAIList *ail, int cLen) {
   //New continueous memory?
    int cLen1=cLen/2, j1, nr, minL = GMAX(64, cLen);
    cLen += cLen1;
    int lenT, len, iter, i, j, k, k0, t;
	for(i=0; i<ail->nctg; i++){
		//1. Decomposition
		TAILstCtg *p    = &ail->ctg[i];
		gdata_t *L1 = p->glist;							//L1: to be rebuilt
		nr 			= p->nr;
		radix_sort_intv(L1, L1+nr);
        if(nr<=minL){
            p->nc = 1, p->lenC[0] = nr, p->idxC[0] = 0;
        }
        else{
        	gdata_t *L0 = malloc(nr*sizeof(gdata_t)); 	//L0: serve as input list
            gdata_t *L2 = malloc(nr*sizeof(gdata_t));   //L2: extracted list
            memcpy(L0, L1, nr*sizeof(gdata_t));
            iter = 0;	k = 0;	k0 = 0;
            lenT = nr;
            while(iter<MAXC && lenT>minL){
                len = 0;
                for(t=0; t<lenT-cLen; t++){
                    uint32_t tt = L0[t].end;
                    j=1;    j1=1;
                    while(j<cLen && j1<cLen1){
                        if(L0[j+t].end>=tt) j1++;
                        j++;
                    }
                    if(j1<cLen1) memcpy(&L2[len++], &L0[t], sizeof(gdata_t));
                    else memcpy(&L1[k++], &L0[t], sizeof(gdata_t));
                }
                memcpy(&L1[k], &L0[lenT-cLen], cLen*sizeof(gdata_t));
                k += cLen, lenT = len;
                p->idxC[iter] = k0;
                p->lenC[iter] = k-k0;
                k0 = k, iter++;
                if(lenT<=minL || iter==MAXC-2){			//exit: add L2 to the end
                    if(lenT>0){
                        memcpy(&L1[k], L2, lenT*sizeof(gdata_t));
                        p->idxC[iter] = k;
                        p->lenC[iter] = lenT;
                        iter++;
                        lenT = 0;	//exit!
                    }
                   	p->nc = iter;
                }
                else memcpy(L0, L2, lenT*sizeof(gdata_t));
            }
            free(L2),free(L0);
        }
        //2. Augmentation
        p->maxE = malloc(nr*sizeof(uint32_t));
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
*/
int32_t TAIList::get_ctg(const char *chr) {
	khint_t k;
	strhash_t *h = (strhash_t*)hc;
	k = kh_get(str, h, chr);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

uint32_t ailist_query(TAIList& ail, int32_t gid, uint32_t qs, uint32_t qe, GVec<int>& r) {
    uint32_t nr = 0; //, m = *mr; //*r = *ir;
    //int32_t gid = get_ctg(ail, chr);
    if(gid>=ail.ctg.Count() || gid<0)return 0;
    r.Clear();
    TAILstCtg *p = ail.ctg[gid];
    for(int k=0; k<p->nc; k++){					//search each component
        int32_t cs = p->idxC[k];
        int32_t ce = cs + p->lenC[k];
        int t;
        if(p->lenC[k]>15){
            t = bSearch(p->glist, cs, ce, qe); 	//rs<qe: inline not better
            if(t>=cs){
		        /*if(nr+t-cs>=m){
		        	m = nr+t-cs + 1024;
		        	r = realloc(r, m*sizeof(uint32_t));
		        }*/
		        while(t>=cs && p->maxE[t]>qs){
		            if(p->glist[t]->end>qs) {
		                r.Add(t);
		                nr++;
		            }
		            t--;
		        }
            }
        }
        else{
        	/*if(nr+ce-cs>=m){
        		m = nr+ce-cs + 1024;
        		r = realloc(r, m*sizeof(uint32_t));
        	}*/
            for(t=cs; t<ce; t++)
                if(p->glist[t]->start<qe && p->glist[t]->end>qs) {
                    //r[nr++] = t;
                	r.Add(t);
                	nr++;
                }
        }
    }
    //*ir = r, *mr = m;
    return nr;
}

