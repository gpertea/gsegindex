CXX   := $(if $(CXX),$(CXX),g++)
GCLIB := ../gclib
INC := -I. -I$(GCLIB)
CXXFLAGS := -g -DNDEBUG -Wall -O3 $(INC) -std=c++11 -fpermissive
LIBS := -lz

LINKER  := $(if $(LINKER),$(LINKER),g++)
DFLAGS := $(if $(LDFLAGS),$(LDFLAGS),-g)

#LIB = AIList.o ailist_main.o
#OBJS = $(addprefix $(OBJ)/, $(LIB))

ifneq (,$(findstring mingw,$(shell ${CXX} -dumpmachine)))
 WINDOWS=1
endif

# File endings
ifdef WINDOWS
 LIBS += -lregex -lws2_32
#else
endif

PROG = ailist gailist
# ailist bedcov bedcov-cpp

all:$(PROG)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

bedcov: bedcov.o cgranges.o $(GCLIB)/GBase.o
		${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
bedcov-cpp: bedcov-iitree.o $(GCLIB)/GBase.o
		${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
bedcov-bfs: bedcov-iitree-bfs.o $(GCLIB)/GBase.o
		${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}

gbedcov: gbedcov.o iutil.o $(GCLIB)/GBase.o $(GCLIB)/GResUsage.o $(GCLIB)/GStr.o
		${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}

#ailist: AIList.co ailist_main.co $(GCLIB)/GBase.o $(GCLIB)/GResUsage.o
ailist: AIList.o iutil.o ailist_main.o $(GCLIB)/GBase.o $(GCLIB)/GResUsage.o
		${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}

gailist: iutil.o gailist_main.o $(GCLIB)/GBase.o $(GCLIB)/GResUsage.o
		${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}

gbedcov.o: iutil.h iutil.cpp gbedcov.cpp
GAIList.o: iutil.h iutil.cpp GAIList.h GRadixSorter.hh
ailist_main.o: iutil.h iutil.cpp AIList.h AIList.cpp
gailist_main.o: iutil.h iutil.cpp GAIList.hh

bedcov-iitree-bfs.o: bedcov-iitree-bfs.cpp IITreeBFS.h
bedcov-iitree.o: bedcov-iitree.cpp IITree.h
AIList.o: AIList.cpp AIList.h

clean:
		rm -fr $(PROG) *.exe *.o a.out *.dSYM $(GCLIB)/GBase.o
