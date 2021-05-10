SRCS = base.cpp msa.cpp model.cpp checkDensity.cpp modelSeq.cpp modelIO.cpp \
	modelGauge.cpp InvertMatrix.cpp PartFunc.cpp ModSeq.cpp MCseq.cpp OptParam.cpp \
	AnaTraj.cpp

OBJS = $(SRCS:.cpp=.o)

BOOST_INC=/opt/local/include
BOOST_LIB=/opt/local/lib
BOOST_LIBPO=boost_program_options-mt
LIBCBLAS=cblas

# BOOST_INC=/home/akinjo/boost/include
# BOOST_LIB=/home/akinjo/boost/lib
# BOOST_LIBPO=boost_program_options
# LIBCBLAS=gslcblas

CC=gcc
CXX=g++
CFLAGS = -O3 -I $(BOOST_INC)
# -Wc++11-extensions 
LIBS= -L $(BOOST_LIB) -l$(BOOST_LIBPO) -lgsl -l$(LIBCBLAS) -lm
.SUFFIXES: .cpp .c .o .hpp

all: lgm_train lgm_trainJ lgm_trainPLM lgm_trainT lgm_mcsamp lgm_wldos lgm_wlsamp lgm_ana lgm_anaJ lgm_pca lgm_hmmemit test
test: test.cpp $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $< $(LIBS)
lgm_train: lgm_train.o $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $< $(LIBS)
lgm_trainPLM: lgm_trainPLM.o $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $< $(LIBS)
lgm_trainT: lgm_trainT.o $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $< $(LIBS)
lgm_trainJ: lgm_trainJ.o $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $< $(LIBS)
lgm_mcsamp: lgm_mcsamp.o $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $< $(LIBS)
lgm_wldos: lgm_wldos.o $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $< $(LIBS)
lgm_wlsamp: lgm_wlsamp.o $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $< $(LIBS)
lgm_ana: lgm_ana.o $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $< $(LIBS)
lgm_hmmemit: lgm_hmmemit.o $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $< $(LIBS)
lgm_anaJ: lgm_anaJ.o $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $< $(LIBS)
lgm_pca: lgm_pca.o $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $< $(LIBS)

.c.o:	
	$(CC) $(CFLAGS) -c $<
.cpp.o:
	$(CXX) $(CFLAGS) -c $<
depend:
	$(CXX) -MM -MG *.cpp > .depend
clean:
	rm -f $(OBJS)
-include .depend
