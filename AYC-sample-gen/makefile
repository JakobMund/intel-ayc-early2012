# using gcc :
COMPILER ?= $(GCC_PATH)g++

# using gcc :
FLAGS ?= -std=c++0x -O3 -Wall $(GCC_SUPPFLAGS)

LDFLAGS ?= -g
EXECUTABLE = ayc-sample-gen

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
	