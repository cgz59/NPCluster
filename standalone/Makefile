CXX = clang++
CXX_FLAGS = -std=c++11 -g -Wall -O2 -DNDEBUG

#R_HOME = `R RHOME` # How to execute?
R_HOME = /Library/Frameworks/R.framework/Resources

R_LIB = -F${R_HOME}/../.. -framework R -Wl,-framework -Wl,CoreFoundation
R_INC_PATH = -I../src -I${R_HOME}/include -I${R_HOME}/library/Rcpp/include
R_LIB_PATH = -L${R_HOME}/lib 

SOURCES = benchmark.cpp
	
OBJS = benchmark.o

%.o:	%.cpp
	$(CXX) $(CXX_FLAGS) -c $*.cpp $(INC) $(R_INC_PATH)
	
benchmark:	$(OBJS)
	$(CXX) $(CXX_FLAGS) -o benchmark $(OBJS) $(R_LIB_PATH) $(R_LIB)


#CXX=ccache clang++


clean:
	/bin/rm -f *.o benchmark


