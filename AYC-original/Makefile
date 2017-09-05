# sample Makefile.
# It compiles every .cpp files in the src/ directory to object files in obj/ directory, and build the ./run executable.
# It automatically handles include dependencies.

# You can modify it as you want (add libraries, compiler flags...).
# But if you want it to properly run on our benchmark service, don't rename or delete variables.

# using icc :
#COMPILER ?= $(ICC_PATH)icpc
# using gcc :
COMPILER ?= $(GCC_PATH)g++

# using icc :
#FLAGS ?= -std=c++0x -U__GXX_EXPERIMENTAL_COMPILER0X__ -xHOST -fast -w1 $(ICC_SUPPFLAGS)
# using gcc :
FLAGS ?= -std=c++0x -O3 -Wall $(GCC_SUPPFLAGS)

LDFLAGS ?= -g
LDLIBS = 
#example if using Intel® Threading Building Blocks :
#LDLIBS = -ltbb -ltbbmalloc

EXECUTABLE = run

SRCS=$(wildcard src/*.cpp)
OBJS=$(SRCS:src/%.cpp=obj/%.o)

all: release

release: $(OBJS)
	$(COMPILER) $(LDFLAGS) -o $(EXECUTABLE) $(OBJS) $(LDLIBS) 

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
	
