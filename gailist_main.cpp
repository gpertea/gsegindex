//=============================================================================
//Read .BED datasets, and then find all overlaps
//by Jianglin Feng  09/05/2018
//Decomposition & simplication: 11/26/2018
//Radix sorting and one-pass loading based on lh3's cgranges: 6/20/2019
//-----------------------------------------------------------------------------
#include "GAIList.hh"
#include "iutil.h"
#define PROGRAM_NAME  "gailist"

int ailist_help(int argc, char **argv, int exit_code);

void loadFromBED(GAIList<int32_t>& a, const char* fn) {
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
			if (ctg) a.add(ctg, st, en, k++);

		}
	} else GError("Error: failed to open file %s\n", fn);
	gzclose(fp);

}

void showHits(GAIList<int32_t>& a) {
   int nh=a.hitCount();
   for (int i=0;i<nh;i++) {
	   const char* ctg=a.hit_ctg(i);
	   uint32_t hstart=a.hit_start(i);
	   uint32_t hend=a.hit_end(i);
	   int32_t d=a.hit(i);
	   printf("hit %d: %s:%d-%d\n", d, ctg, hstart, hend);
   }
}

int main(int argc, char **argv)
{
    int cLen = 20, pmode = 0;//print mode: 0 print hitd[i]; 1: print total; 2: print components
    for(int i=3; i<argc; i++){
    	if(strcmp(argv[i], "-L")==0 && i+1<argc)
    		cLen = atoi(argv[i+1]);
    	else if(strcmp(argv[i], "-P")==0 && i+1<argc)
    		pmode = atoi(argv[i+1]);
    }
	if(argc<3)
        return ailist_help(argc, argv, 0);

   	//clock_t start, end1, end2, end3;
    //start = clock();

    //1. Read interval data
  GResUsage ru;
  //ru.start();

  //AIList *ail =  greadBED(argv[1]);
  GAIList<int32_t> gail;
  loadFromBED(gail, argv[1]);

   //end1 = clock();
    //printf("loading time: %f\n", ((double)(end1-start))/CLOCKS_PER_SEC);
    //2. Construct ailist

   //gailist_construct(ail, cLen);
   gail.build(cLen);

   //ru.stop();
   //double mtime=ru.elapsed()/1000000; // in seconds
   //double memused=ru.memoryUsed();
   //GMessage("%s loaded (in %.2f sec, using %.2f KB)\n", argv[1], mtime, memused);
   //ru.start();
   /*
   if(pmode==2){
  	  for(uint i=0;i<ail.ctgs.Count();i++){
			 AICtgData& p = ail.ctgs[i];
			 printf("%s\tnr= %lld, nc=%i\n", p.ctg, (long long)p.nr, p.nc);
			 for(uint j=0;j<p.nc;j++)
				  printf(" %i\t%i  \n", p.idxC[j], p.lenC[j]);
		  }
    }
    */
    //end2 = clock();
    //printf("constru time: %f\n", ((double)(end2-end1))/CLOCKS_PER_SEC);
    //ru.stop();
    //mtime=ru.elapsed()/1000000; // in seconds
    //GMessage("%s AIList built in %.2f sec\n", argv[1], mtime, memused);
    //3. Search
	int64_t nol = 0;
	uint32_t nhits=0;//, mr=1000000; //mr: initial capacity of hits[]
	gzFile fp;
	if ((fp=gzopen(argv[2], "r"))) {
		GFStream<gzFile, int (*)(gzFile, voidp, unsigned int)> fs(fp, gzread);
		Gcstr line;
		while (fs.getUntil(fs.SEP_LINE, line)>=0) {
			if (line.len()==0) continue;
			uint32_t st1, en1;
			char *ctg;
			ctg = parse_bed(line(), &st1, &en1);
			if (ctg == 0) continue;
			//nhits = gailist_query(ail, ctg, st1, en1, &mr, &hits);
			nhits=gail.query(ctg, st1, en1);
			//nhits=ail.query(ctg, st1, en1, hits);
			if (pmode==0 && nhits>0)
				printf("%s\t%d\t%d\t%d\n", ctg, st1, en1, nhits);
			//showHits(gail);
			nol += nhits;
		}
	} else GError("Error: failed to open file %s\n", argv[2]);

	gzclose(fp);

    return 0;
}

int ailist_help(int argc, char **argv, int exit_code)
{
    fprintf(stderr,"gailist\n" "usage:   gailist database-file(.bed) query-file(.bed) [-L coverage-length] \n");
    return exit_code;
}

