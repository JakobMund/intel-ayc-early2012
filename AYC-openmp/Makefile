# using icc :
#COMPILER ?= $(ICC_PATH)icpc
# using gcc :
COMPILER ?= $(GCC_PATH)g++

# using icc :
# FLAGS ?= -std=c++0x -U__GXX_EXPERIMENTAL_COMPILER0X__ -xHOST -fast -w1 $(ICC_SUPPFLAGS)
# advanced flags using (icc):
#  -xSSE4.2: enable SSE4.2 instructions 
#FLAGS ?= -std=c++0x -U__GXX_EXPERIMENTAL_COMPILER0X__ -xHOST -fast -w1 -pthread $(ICC_SUPPFLAGS)

# using gcc :
#FLAGS ?= -std=c++0x -O3 -Wall $(GCC_SUPPFLAGS)
# advanced flags using (gcc):
#  -0fast: enables even more aggressive optimizations, even those not compliying to the standard (this includes -ffast-math)
#  -march=corei7-avs: optimize for the core i3/i5/i7 machine architecture; should include the xeon processors
#  -march=native & -mtune=native: uses and optimizes all instruction available on the machine this program is compiled on
#  -fprefetch-loop-arrays: optimize loops on large arrays by trying to prefetch data
#  -msse4.2: enable instructions offered by SSE4.2
#FLAGS ?= -std=c++0x -Ofast -Wall -march=corei7 -mtune=corei7 $(GCC_SUPPFLAGS)
FLAGS ?= -fopenmp -std=c++0x -O3 -Wall -march=native -mtune=native -ffast-math -fprefetch-loop-arrays -msse4.2 $(GCC_SUPPFLAGS)
FLAGSA ?= -c -g -Wa,-a,-ad

LDFLAGS ?= -g -fopenmp
EXECUTABLE = run

SRCS=$(wildcard src/*.cpp)
OBJS=$(SRCS:src/%.cpp=obj/%.o)

all: release

release: $(OBJS)
	$(COMPILER) $(LDFLAGS) -o $(EXECUTABLE) $(OBJS) 

obj/%.o: src/%.cpp
	$(COMPILER) $(FLAGS) -o $@ -c $<

clean:
	rm -f obj/*

dist-clean: clean
	rm -f $(EXECUTABLE) *~ .depend
	
#automatically handle include dependencies
depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	@$(foreach SRC, $(SRCS), $(COMPILER) $(FLAGS) -MT $(SRC:src/%.cpp=obj/%.o) -MM $(SRC) >> .depend;)

include .depend
	
