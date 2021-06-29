CXX   := $(if $(CXX),$(CXX),g++)
GCLIB:=../gclib
INC:= -I. -I$(GCLIB)
CXXFLAGS= -g -Wall -O2 $(INC) -std=c++11 -fpermissive
LIBS=-lz

LINKER  := $(if $(LINKER),$(LINKER),g++)
DFLAGS := $(if $(LDFLAGS),$(LDFLAGS),-g)

#LIB = AIList.o ailist_main.o
#OBJS = $(addprefix $(OBJ)/, $(LIB))

EXE=bedcov bedcov-bfs ailist


all:$(EXE)

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
		rm -fr *.o a.out *.dSYM $(GCLIB)/GBase.o $(EXE)
