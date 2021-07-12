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

AIList *gailist_init(void)
{
	AIList *ail = malloc(1*sizeof(AIList));
	//ail->hc = kh_init(str);
	ail->ctghash=new GHashMap<const char*, int32_t>();
	ail->ctghash->resize(64);
	ail->nctg = 0;
	ail->mctg = 32;
	ail->ctg = malloc(ail->mctg*sizeof(ctg_t));
	return ail;
}

void GAIList::init() {
	hc = kh_init(str);
	nctg = 0;
	mctg = 32;
	ctg = malloc(mctg*sizeof(ctg_t));
}


void gailist_destroy(AIList *ail)
{
	int32_t i;
	if (ail == 0) return;
	for (i = 0; i < ail->nctg; ++i){
		free(ail->ctg[i].name);
		free(ail->ctg[i].glist);
		free(ail->ctg[i].maxE);
	}
	free(ail->ctg);
	//kh_destroy(str, (strhash_t*)ail->hc);
	delete ail->ctghash;
	free(ail);
}

void GAIList::destroy() {
	for (int i = 0; i < nctg; ++i){
		free(ctg[i].name);
		free(ctg[i].glist);
		free(ctg[i].maxE);
	}
	free(ctg);
	kh_destroy(str, (strhash_t*)hc);
}



void gailist_add(AIList *ail, const char *chr, uint32_t s, uint32_t e, uint32_t v)
{
	if(s > e) return;
	bool cnew=false;
	uint64_t hidx=ail->ctghash->addIfNew(chr, ail->nctg, cnew);
	ctg_t *q;
	if (cnew) { //new contig
		if (ail->nctg == ail->mctg)
			{ EXPAND(ail->ctg, ail->mctg); }
		q = &ail->ctg[ail->nctg++];
		q->name=strdup(chr);
		q->nr=0; q->mr=64;
		GMALLOC( q->glist, (q->mr * sizeof(AIData)) );
		ail->ctghash->setKey(hidx, q->name);
	} else {
		q = &ail->ctg[ail->ctghash->getValue(hidx)];
	}

	if (q->nr == q->mr)
		{ EXPAND(q->glist, q->mr); }
	AIData *p = &q->glist[q->nr++];
	p->start = s;
	p->end   = e;
	p->didx = v;
}

void GAIList::add(const char *chr, uint32_t s, uint32_t e, uint32_t payload) {
	if(s > e) return;
	bool cnew=false;
	uint64_t hidx=ctghash->addIfNew(chr, nctg, cnew);
	ctg_t *q;
	if (cnew) { //new contig
		if (ail->nctg == ail->mctg)
			{ EXPAND(ail->ctg, ail->mctg); }
		q = &ail->ctg[ail->nctg++];
		q->name=strdup(chr);
		q->nr=0; q->mr=64;
		GMALLOC( q->glist, (q->mr * sizeof(AIData)) );
		ail->ctghash->setKey(hidx, q->name);
	} else {
		q = &ail->ctg[ail->ctghash->getValue(hidx)];
	}

	if (q->nr == q->mr)
		{ EXPAND(q->glist, q->mr); }
	AIData *p = &q->glist[q->nr++];
	p->start = s;
	p->end   = e;
	p->didx = v;
}

//-------------------------------------------------------------------------------
AIList* greadBED(const char* fn)
{   //faster than strtok()
	gzFile fp;
	AIList *ail;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	int32_t k = 0;
	if ((fp = gzopen(fn, "r")) == 0)
		return 0;
	ks = ks_init(fp);
	ail = gailist_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		char *ctg;
		int32_t st, en;
		ctg = parse_bed(str.s, &st, &en);
		if (ctg) gailist_add(ail, ctg, st, en, k++);
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	return ail;
}

void gailist_construct(AIList *ail, int cLen)
{
    int cLen1=cLen/2, j1, nr, minL = GMAX(64, cLen);
    cLen += cLen1;
    int lenT, len, iter, i, j, k, k0, t;
	for(i=0; i<ail->nctg; i++){
		//1. Decomposition
		ctg_t *p    = &ail->ctg[i];
		AIData *L1 = p->glist;							//L1: to be rebuilt
		nr 			= p->nr;
		radix_sort_intv(L1, L1+nr);
        if(nr<=minL){
            p->nc = 1, p->lenC[0] = nr, p->idxC[0] = 0;
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
                k += cLen, lenT = len;
                p->idxC[iter] = k0;
                p->lenC[iter] = k-k0;
                k0 = k, iter++;
                if(lenT<=minL || iter==MAXC-2){			//exit: add L2 to the end
                    if(lenT>0){
                        memcpy(&L1[k], L2, lenT*sizeof(AIData));
                        p->idxC[iter] = k;
                        p->lenC[iter] = lenT;
                        iter++;
                        lenT = 0;						//exit!
                    }
                   	p->nc = iter;
                }
                else memcpy(L0, L2, lenT*sizeof(AIData));
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
/*
//-------------------------------------------------------------------------------
GAIList::GAIList(uint max):ctgs(max) { hc = kh_init(str); } //init
GAIList::~GAIList() { //ailist_destroy
    	//ctglst is going to clear itself
	    //but ctgs won't
	   for (int i=0;i<ctgs.Count();++i)
		   ctgs[i].destroy();
    	kh_destroy(str, (strhash_t*)this->hc);
}
*/
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

/*
void GAIList::add(const char *chr, uint32_t s, uint32_t e, uint32_t payload) {
	if(s > e)return;
	int absent;
	khint_t k;
	strhash_t *h = (strhash_t*)(this->hc);
	k = kh_put(str, h, chr, &absent);
	AICtgData* ctg;
	if (absent) { //add new AICtgData
		//kh_val(h, k) = ctglst.Count();
		//AICtgData* newctg=new AICtgData(chr);
		ctg=&(ctgs.Add());
		ctg->init(chr);
		kh_val(h, k) = ctgs.Count()-1;
		kh_key(h, k) = ctg->ctg;
	} else {
		ctg = & ctgs[kh_val(h, k)];
	}
	//q.reglist.Add({s, e, payload});
	AIData& p = ctg->add();
	p.start = s; p.end = e; p.didx = payload;
}
*/

// TFunc for GFStream:
//typedef int (*ft_gzread)(gzFile, voidp, unsigned int);
/*
AIList* readBED(const char* fn)
{   //faster than strtok()
	gzFile fp;
	AIList *ail;
	//int32_t k = 0;
	ail = ailist_init();
	if ((fp=gzopen(fn, "r"))) {
		GFStream<gzFile, int (*)(gzFile, voidp, unsigned int)> fs(fp, gzread);
		Gcstr line;
		while (fs.getUntil(fs.SEP_LINE, line)>=0) {
			if (line.len()==0) continue;
			char *ctg;
			int32_t stailist, en;
			ctg = parse_bed(line(), &st, &en);
			if (ctg) ailist_add(ail, ctg, st, en);//, k++);

		}
	} else GError("Error: failed to open file %s\n", fn);
	gzclose(fp);
	return ail;
}
*/
/*
int32_t get_ctg(const AIList* ail, const char *chr)
{
	khint_t k;
	strhash_t *h = (strhash_t*)ail->hc;
	k = kh_get(str, h, chr);
	return k == kh_end(h)? -1 : kh_val(h, k);
}
*/

int32_t get_ctg(const AIList* ail, const char *chr) {
	int32_t* hp=ail->ctghash->Find(chr);
	return hp ? *hp : -1;
}

uint32_t gailist_query(AIList *ail, char *chr, uint32_t qs, uint32_t qe, uint32_t *mr, uint32_t **ir)
{
    uint32_t nr = 0, m = *mr, *r = *ir;
    int32_t gid = get_ctg(ail, chr);
    if(gid>=ail->nctg || gid<0)return 0;
    ctg_t *p = &ail->ctg[gid];
    for(int k=0; k<p->nc; k++){					//search each component
        int32_t cs = p->idxC[k];
        int32_t ce = cs + p->lenC[k];
        int32_t t;
        if(p->lenC[k]>15){
            t = bSearch(p->glist, cs, ce, qe); 	//rs<qe: inline not better
            if(t>=cs){
		        if(nr+t-cs>=m){
		        	m = nr+t-cs + 1024;
		        	r = realloc(r, m*sizeof(uint32_t));
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
        		r = realloc(r, m*sizeof(uint32_t));
        	}
            for(t=cs; t<ce; t++)
                if(p->glist[t].start<qe && p->glist[t].end>qs)
                    r[nr++] = t;
        }
    }
    *ir = r, *mr = m;
    return nr;
}
/*
void GAIList::build(int cLen) {
	int cLen1=cLen/2;
	uint32_t minL = GMAX(64, cLen);
	uint32_t nr;
	cLen += cLen1;
	uint32_t lenT, len, iter, j, t, k, k0;
	for(uint i=0; i<ctgs.Count(); i++){ //for each chromosome
		//1. Decomposition
		AICtgData& ctgdta = ctgs[i];
		AIData* L1 = ctgdta.reglist; //L1: to be rebuilt
		nr = ctgdta.nr;
		radix_sort_intv(L1, L1+nr);
		if (nr<=minL) {
			ctgdta.nc = 1;
			ctgdta.lenC[0] = nr;
			ctgdta.idxC[0] = 0;
		}
		else{
			AIData *L0 = (AIData*) malloc(nr*sizeof(AIData)); 	//L0: serve as input list
			AIData *L2 = (AIData*) malloc(nr*sizeof(AIData));   //L2: extracted list
			//----------------------------------------
			AIData *D0 = (AIData*) malloc(nr*sizeof(AIData)); 	//D0:
			int32_t *di = (int32_t*) malloc(nr*sizeof(int32_t));	//int64_t?
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
*/

int32_t GAIList::getCtg(const char *chr) {
	strhash_t *h = (strhash_t*)(this->hc);
	khint_t k = kh_get(str, h, chr);
	return k == kh_end(h) ? -1 : kh_val(h, k);
}

/*
uint32_t GAIList::query(char *chr, uint32_t qs, uint32_t qe, GDynArray<uint32_t>& hits) {
    //uint32_t nr = 0, m = *mr, *r = *ir;
    //int32_t gid = this->getCtg(chr);
	int32_t gid = get_ctg(this->hc, chr);

    if (gid<0 || (uint32_t)gid >= ctgs.Count() )
    	return 0; //no such contig
    AICtgData& p = ctgs[gid];
    for(uint k=0; k<p.nc; k++){     //search each component
        int32_t cs = p.idxC[k];
        int32_t ce = cs + p.lenC[k];
        int32_t t;
        if(p.lenC[k]>15) { // binary search
            t = bSearch(p.reglist, cs, ce, qe);  //rs<qe: inline not better
            if(t>=cs){
		        while(t>=cs && p.maxE[t]>qs){
		            if(p.reglist[t].end>qs)
		                hits.Add(t);
		            t--;
		        }
            }
        }
        else{ //linear search
            for(t=cs; t<ce; t++)
                if(p.reglist[t].start<qe && p.reglist[t].end>qs)
                	hits.Add(t);
        }
    }
    return hits.Count();
}
*/
