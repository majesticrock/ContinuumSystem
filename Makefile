CXX = mpicxx

# -isystem rather than -I as this supresses warnings that occur during
# the compilation of the eigen library (I cant fix them anyways)
INCLUDEFLAGS = -isystem ~/usr/local/include

CXXFLAGS = $(WARNINGS) -std=c++20 $(OPT) -fopenmp

LDLIBS = -L/sw/gcc/5.3.0/rtf/lib64 -L/home/joshua/usr/local/include/boost_lib/ -lboost_serialization -lboost_iostreams -lz

WARNINGS = -Wall -Wno-sign-compare

OPT = -march=native -O3# -ffast-math

#COMMUTE_SRCS=Coefficient.cpp IndexWrapper.cpp Momentum.cpp Operator.cpp Term.cpp WickCleaner.cpp WickOperator.cpp WickOperatorTemplate.cpp WickTerm.cpp 
UTIL_SRCS=Utility/InputFileReader.cpp

CONT_SRCS=SCModel.cpp
PART_SRCS=ContinuumSystem.cpp
SRCS=$(addprefix SymbolicOperators/, $(COMMUTE_SRCS)) $(addprefix Continuum/, $(CONT_SRCS)) $(PART_SRCS) $(UTIL_SRCS)

OBJS=$(addprefix build/, $(subst .cpp,.o,$(SRCS)))

all: sources/SymbolicOperators sources/Utility build build/main 

debug: CXXFLAGS += -g -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment
debug: build build/main

build/main: $(OBJS) | build 
	$(CXX) $(INCLUDEFLAGS) -o build/main $(OBJS) $(CXXFLAGS) $(LDLIBS)

build/%.o: sources/%.cpp# sources/%.hpp
	$(CXX) $(INCLUDEFLAGS) $< -o $@ -c $(CXXFLAGS)

sources/Utility:
	ln -s ../../Utility/sources sources/Utility

sources/SymbolicOperators:
	ln -s ../../FermionCommute/sources sources/SymbolicOperators

build:
	mkdir -p build
	mkdir -p build/Continuum
	mkdir -p build/Utility
	mkdir -p build/SymbolicOperators

clean:
	rm -rf build
