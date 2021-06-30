CXX   := $(if $(CXX),$(CXX),g++)
GCLIB := ../gclib
INC := -I. -I$(GCLIB)
CXXFLAGS := -g -Wall -O2 $(INC) -std=c++11 -fpermissive
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

PROG=bedcov bedcov-bfs ailist

all:$(PROG)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

bedcov: bedcov-iitree.o $(GCLIB)/GBase.o
		${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
bedcov-bfs: bedcov-iitree-bfs.o $(GCLIB)/GBase.o
		${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}

ailist: AIList.o ailist_main.o $(GCLIB)/GBase.o
		${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}


bedcov-iitree-bfs.o: bedcov-iitree-bfs.cpp IITreeBFS.h
bedcov-iitree.o: bedcov-iitree.cpp IITree.h
AIList.o: AIList.cpp AIList.h


clean:
		rm -fr $(PROG) *.exe *.o a.out *.dSYM $(GCLIB)/GBase.o
