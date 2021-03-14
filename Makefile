CXX   := $(if $(CXX),$(CXX),g++)
CXXFLAGS= -g -Wall -std=c++11 -fpermissive
CFLAGS=-g -Wall -Wc++-compat
LIBS=-lz

LINKER  := $(if $(LINKER),$(LINKER),g++)
DFLAGS := $(if $(LDFLAGS),$(LDFLAGS),-g)

#LIB = AIList.o ailist_main.o
#OBJS = $(addprefix $(OBJ)/, $(LIB))

EXE=bedcov ailist


all:$(EXE)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

bedcov: bedcov-iitree.o
		${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
ailist: AIList.o ailist_main.o
		${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}


bedcov-iitree.o: bedcov-iitree.cpp IITreeBFS.h
AIList.o: AIList.cpp AIList.h


clean:
		rm -fr *.o a.out *.dSYM $(EXE)
