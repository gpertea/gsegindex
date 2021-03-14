#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "IITreeBFS.h"
#include "khash.h"
#include "kseq.h"

#include "GBase.h"
#include "GVec.hh"

KSTREAM_INIT(gzFile, gzread, 0x10000)

typedef IITree<int32_t, int32_t> ITree;

KHASH_MAP_INIT_STR(idx, ITree*)

char *parse_bed(char *s, int32_t *st_, int32_t *en_)
{
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

khash_t(idx) *read_bed(const char *fn)
{
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	khash_t(idx) *h;

	if ((fp = gzopen(fn, "r")) == 0)
		return 0;
	ks = ks_init(fp);
	h = kh_init(idx);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		char *ctg;
		int32_t st, en;
		ctg = parse_bed(str.s, &st, &en);
		if (ctg) {
			khint_t k;
			int absent;
			k = kh_put(idx, h, ctg, &absent);
			if (absent) {
				kh_key(h, k) = strdup(ctg);
				kh_val(h, k) = new ITree;
			}
			kh_val(h, k)->add(st, en, kh_val(h, k)->size());
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);

	for (khint_t k = 0; k < kh_end(h); ++k)
		if (kh_exist(h, k))
			kh_val(h, k)->index();
	return h;
}

int main(int argc, char *argv[])
{
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};

	if (argc < 3) {
		printf("Usage: bedcov <loaded.bed> <streamed.bed>\n");
		return 0;
	}

	khash_t(idx) *h = read_bed(argv[1]);

	fp = gzopen(argv[2], "r");
	assert(fp);
	ks = ks_init(fp);
	std::vector<size_t> a; //it will receive the results
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		int32_t st1, en1;
		ITree *tree=NULL;
		char *ctg;
		int numhits=0;
		ctg = parse_bed(str.s, &st1, &en1);
		if (ctg == 0) continue;
		a.clear();
		khint_t k = kh_get(idx, h, ctg);
		if (k != kh_end(h)) {
		  tree = kh_val(h, k);
		  tree->overlap(st1, en1, a);
		}
		printf(">Qry_%s:[%d-%d] (%d hits)\n", ctg, st1, en1, (int)a.size());
		if (a.size()==0) continue;
		int32_t cnt = 0, cov = 0, cov_st = 0, cov_en = 0;
		GVec<GSeg> rlst(a.size());
		for (size_t j = 0; j < a.size(); ++j) {
			//int32_t st0 = tree->start(a[j]), en0 = tree->end(a[j]);
			GSeg seg(tree->start(a[j]), tree->end(a[j]));
			rlst.Add(seg);
			continue;
			/*
			if (st0 < st1) st0 = st1;
			if (en0 > en1) en0 = en1;
			if (st0 > cov_en) {
				cov += cov_en - cov_st;
				cov_st = st0, cov_en = en0;
			} else cov_en = cov_en > en0? cov_en : en0;
			++cnt;
			*/
		}
		rlst.Sort();
		for (int i=0;i<rlst.Count();++i)
		  printf("%s\t%d\t%d\n", ctg, rlst[i].start, rlst[i].end);
		//cov += cov_en - cov_st;
		//printf("%s\t%d\t%d\t%d\t%d\n", ctg, st1, en1, cnt, cov);
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);

	for (khint_t k = 0; k < kh_end(h); ++k)
		if (kh_exist(h, k)) {
			delete kh_val(h, k);
			free((void*)kh_key(h, k));
		}
	kh_destroy(idx, h);
	return 0;
}
