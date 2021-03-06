//=============================================================================
//Read .BED datasets, and then find all overlaps
//by Jianglin Feng  09/05/2018
//Decomposition & simplication: 11/26/2018
//Radix sorting and one-pass loading based on lh3's cgranges: 6/20/2019
//-----------------------------------------------------------------------------
#include "AIList.h"
#define PROGRAM_NAME  "ailist"
#define MAJOR_VERSION "0"
#define MINOR_VERSION "1"
#define REVISION_VERSION "1"
#define BUILD_VERSION "0"
#define VERSION MAJOR_VERSION "." MINOR_VERSION "." REVISION_VERSION

#include "GVec.hh"

int ailist_help(int argc, char **argv, int exit_code) {
    fprintf(stderr,"%s, v%s\n" "usage:   %s database-file(.bed) query-file(.bed)\n",
            PROGRAM_NAME, VERSION, PROGRAM_NAME);
    return exit_code;
}

int main(int argc, char **argv) {
    //int cLen = 20, pmode = 0;//print mode: 0 print hitd[i]; 1: print total; 2: print components
	if(argc<3)
        return ailist_help(argc, argv, 0);
    //1. Read interval data
    TAIList ail;
    readBED(ail, argv[1]);
    //printf("loading time: %f\n", ((double)(end1-start))/CLOCKS_PER_SEC);

    //2. Construct ailist
    ail.build(10);
    //3. Search
	int64_t nol = 0;
	uint32_t nhits=0; //, mr=100;  //, m=0;
	//uint32_t *hits=malloc(mr*sizeof(uint32_t));
    GVec<int> hits;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	gzFile fp = gzopen(argv[2], "r");
	assert(fp);
	ks = ks_init(fp);
	int32_t st0=0, en0=0, st1, en1;
	char *ctg;
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		ctg = parse_bed(str.s, &st1, &en1);
		if (ctg == 0) continue;
		nhits=0;
		int32_t ctgIdx=ail.get_ctg(ctg);
		if (ctgIdx>=0) {
		  nhits = ailist_query(ail, ctgIdx, st1, en1, hits);
		}
		//GMessage(">>> query interval [%d, %d] (%d results)\n", st1, en1, nhits);
		printf(">Qry_%s:[%d-%d] (%d hits)\n", ctg, st1, en1, nhits);
		if (nhits>0) {
			GVec<GSeg> rlst(nhits);
			for (uint h=0;h<nhits;h++) {
			  GSeg &hd = *(ail.ctg[ctgIdx]->glist[hits[h]]);
			  GSeg seg(hd.start, hd.end);
			  rlst.Add(seg);
			  //printf("%s\t%d\t%d\n", ctg, hd.start, hd.end);
			}
			rlst.Sort();
			for (int i=0;i<rlst.Count();++i)
				printf("%s\t%d\t%d\n", ctg, rlst[i].start, rlst[i].end);
		}
		nol += nhits;
		//if(pmode==1 && nhits>0)
		// printf("%i\t%s:\t %i\t %i\t %ld\n", m++, ctg, st1, en1, (long)nhits);
	}
	//end3 = clock();
    //printf("Total %lld\n", (long long)nol);
    //printf("query time: %f\n", ((double)(end3-end2))/CLOCKS_PER_SEC);

  	free(str.s);
	//free(hits);
	gzclose(fp);
	ks_destroy(ks);
	//ailist_destroy(ail);
    return 0;
}

