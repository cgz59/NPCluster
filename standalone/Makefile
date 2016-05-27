CXX = clang++
CXX_FLAGS = -std=c++11 -g -Wall -O2 -DNDEBUG -mavx #-fverbose-asm

#R_HOME = `R RHOME` # How to execute?
R_HOME = /Library/Frameworks/R.framework/Resources
R_LIB_INC = /Users/msuchard/Library/R/3.2
#R_HOME = `Rscript -r "R.home()"`

R_LIB = -F${R_HOME}/../.. -framework R -Wl,-framework -Wl,CoreFoundation -lRInside -ltbb
R_INC_PATH = -I../src -I${R_HOME}/include \
             -I${R_LIB_INC}/library/Rcpp/include \
             -I${R_LIB_INC}/library/BH/include \
             -I${R_LIB_INC}/library/RInside/include \
             -I${R_LIB_INC}/library/RcppNT2/include
             
R_LIB_PATH =	-L${R_LIB_INC}/library/RInside/lib \
				-L${R_HOME}/lib \
				-L${R_LIB_INC}/library/RcppParallel/lib

SOURCES = benchmark.cpp
	
OBJS = benchmark.o

%.o:	%.cpp
	$(CXX) $(CXX_FLAGS) -c $*.cpp $(INC) $(R_INC_PATH)
	
benchmark:	$(OBJS)
	$(CXX) $(CXX_FLAGS) -o benchmark $(OBJS) $(R_LIB_PATH) $(R_LIB)

assembly:	$(OBJS)
	$(CXX) $(CXX_FLAGS) -c -S benchmark.cpp $(INC) $(R_INC_PATH)

clean:
	/bin/rm -f *.o benchmark


