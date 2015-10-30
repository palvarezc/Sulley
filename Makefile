CXX      = g++

# A list of directories
BASE_DIR = $(shell pwd)
LIB_DIR  = $(BASE_DIR)/lib
SRC_DIR  = $(BASE_DIR)/src
OBJ_DIR  = $(BASE_DIR)/obj
INC_DIR  = $(BASE_DIR)/include
BIN_DIR = $(BASE_DIR)/bin

LCGBOOST = /usr/local
BOOST_INC_DIR  = -I$(LCGBOOST)/include/
BOOST_LIB_OPT  = -L$(LCGBOOST)/lib -lboost_filesystem

ROOT_FLAGS   = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOT_LIBS    = $(shell $(ROOTSYS)/bin/root-config --libs) -lTreePlayer -lTreeViewer -lMinuit -lXMLIO -lMLP -lRIO
ROOFIT_LIBS  = ${ROOT_LIBS} -lRooFitCore -lRooFit -lRooStats

FLAGS   = -O3 -Wall -fPIC -D__ROOFIT_NOBANNER -std=c++0x
FLAGS       += ${ROOT_FLAGS} -lRooFitCore -lRooFit -lRooStats ${BOOST_LIB_OPT} ${BOOST_INC_DIR}
FLAGS       += -I$(INC_DIR)

# all: toystudy2DHistBremCatTsallisBkg # mainTestTsallis toystudy1DHistBremCatExpoBkg toystudy1DHistBremCatExpoBkgMHCut mainTestCB toystudy2DHOPHistBremCat 
all: $(BIN_DIR)/maintest $(BIN_DIR)/toystudy $(BIN_DIR)/toystudy2DHistBremCatTsallisBkg $(BIN_DIR)/toystudy1DHistBremCatExpoBkg 

$(BIN_DIR)/maintest : $(OBJ_DIR)/maintest.o $(OBJ_DIR)/RooMcorMvisTsallis.o $(OBJ_DIR)/usefulFunctions.o
	$(CXX) $(FLAGS) -o $(BIN_DIR)/maintest $(OBJ_DIR)/maintest.o $(OBJ_DIR)/RooMcorMvisTsallis.o $(OBJ_DIR)/usefulFunctions.o $(ROOFIT_LIBS)

$(OBJ_DIR)/maintest.o : $(SRC_DIR)/maintest.cc $(INC_DIR)/RooMcorMvisTsallis.h $(INC_DIR)/usefulFunctions.h
	$(CXX) $(FLAGS) -c $(SRC_DIR)/maintest.cc -o $(OBJ_DIR)/maintest.o $(ROOFIT_LIBS)


$(BIN_DIR)/toystudy : $(OBJ_DIR)/toystudy.o $(OBJ_DIR)/RooMcorMvisTsallis.o $(OBJ_DIR)/usefulFunctions.o $(OBJ_DIR)/fitter_utils.o
	$(CXX) $(FLAGS) -o $(BIN_DIR)/toystudy $(OBJ_DIR)/toystudy.o $(OBJ_DIR)/RooMcorMvisTsallis.o $(OBJ_DIR)/usefulFunctions.o $(OBJ_DIR)/fitter_utils.o $(ROOFIT_LIBS)

$(OBJ_DIR)/toystudy.o : $(SRC_DIR)/toystudy.cc $(INC_DIR)/RooMcorMvisTsallis.h $(INC_DIR)/usefulFunctions.h $(INC_DIR)/fitter_utils.h
	$(CXX) $(FLAGS) -c $(SRC_DIR)/toystudy.cc -o $(OBJ_DIR)/toystudy.o $(ROOFIT_LIBS)


$(BIN_DIR)/toystudy2DHistBremCatTsallisBkg : $(OBJ_DIR)/toystudy2DHistBremCatTsallisBkg.o $(OBJ_DIR)/RooMcorMvisTsallis.o $(OBJ_DIR)/usefulFunctions.o $(OBJ_DIR)/fitter_utils.o
	$(CXX) $(FLAGS) -o $(BIN_DIR)/toystudy2DHistBremCatTsallisBkg $(OBJ_DIR)/toystudy2DHistBremCatTsallisBkg.o $(OBJ_DIR)/RooMcorMvisTsallis.o $(OBJ_DIR)/usefulFunctions.o $(OBJ_DIR)/fitter_utils.o $(ROOFIT_LIBS)

$(OBJ_DIR)/toystudy2DHistBremCatTsallisBkg.o : $(SRC_DIR)/toystudy2DHistBremCatTsallisBkg.cc $(INC_DIR)/RooMcorMvisTsallis.h $(INC_DIR)/usefulFunctions.h $(INC_DIR)/fitter_utils.h
	$(CXX) $(FLAGS) -c $(SRC_DIR)/toystudy2DHistBremCatTsallisBkg.cc -o $(OBJ_DIR)/toystudy2DHistBremCatTsallisBkg.o $(ROOFIT_LIBS)


$(BIN_DIR)/toystudy1DHistBremCatExpoBkg : $(OBJ_DIR)/toystudy1DHistBremCatExpoBkg.o $(OBJ_DIR)/usefulFunctions.o 
	$(CXX) $(FLAGS) -o $(BIN_DIR)/toystudy1DHistBremCatExpoBkg $(OBJ_DIR)/toystudy1DHistBremCatExpoBkg.o $(OBJ_DIR)/usefulFunctions.o $(ROOFIT_LIBS)

$(OBJ_DIR)/toystudy1DHistBremCatExpoBkg.o : $(SRC_DIR)/toystudy1DHistBremCatExpoBkg.cc $(INC_DIR)/usefulFunctions.h 
	$(CXX) $(FLAGS) -c $(SRC_DIR)/toystudy1DHistBremCatExpoBkg.cc -o $(OBJ_DIR)/toystudy1DHistBremCatExpoBkg.o $(ROOFIT_LIBS)


# toystudy2DHOPHistBremCat : toystudy2DHOPHistBremCat.o usefulFunctions.o
# 	$(CXX) $(FLAGS) -o toystudy2DHOPHistBremCat toystudy2DHOPHistBremCat.o usefulFunctions.o $(ROOFIT_LIBS)

# toystudy2DHOPHistBremCat.o : toystudy2DHOPHistBremCat.cc usefulFunctions.h
# 	$(CXX) $(FLAGS) -c toystudy2DHOPHistBremCat.cc -o toystudy2DHOPHistBremCat.o $(ROOFIT_LIBS)


# toystudy1DHistBremCatExpoBkgMHCut : toystudy1DHistBremCatExpoBkgMHCut.o usefulFunctions.o
# 	$(CXX) $(FLAGS) -o toystudy1DHistBremCatExpoBkgMHCut toystudy1DHistBremCatExpoBkgMHCut.o usefulFunctions.o $(ROOFIT_LIBS)

# toystudy1DHistBremCatExpoBkgMHCut.o : toystudy1DHistBremCatExpoBkgMHCut.cc usefulFunctions.h
# 	$(CXX) $(FLAGS) -c toystudy1DHistBremCatExpoBkgMHCut.cc -o toystudy1DHistBremCatExpoBkgMHCut.o $(ROOFIT_LIBS)


# $(BIN_DIR)/mainTestTsallis : $(OBJ_DIR)/mainTestTsallis.o $(OBJ_DIR)/RooMcorMvisTsallis.o
# 	$(CXX) $(FLAGS) -o $(BIN_DIR)/mainTestTsallis $(OBJ_DIR)/mainTestTsallis.o $(OBJ_DIR)/RooMcorMvisTsallis.o $(ROOFIT_LIBS)

# $(OBJ_DIR)/mainTestTsallis.o : $(SRC_DIR)/mainTestTsallis.cc $(INC_DIR)/RooMcorMvisTsallis.h
# 	$(CXX) $(FLAGS) -c $(SRC_DIR)/mainTestTsallis.cc -o $(OBJ_DIR)/mainTestTsallis.o $(ROOFIT_LIBS)


# mainTestCB : mainTestCB.o MyCB.o
# 	$(CXX) $(FLAGS) -o mainTestCB mainTestCB.o MyCB.o $(ROOFIT_LIBS)

# mainTestCB.o : mainTestCB.cc MyCB.h
# 	$(CXX) $(FLAGS) -c mainTestCB.cc -o mainTestCB.o $(ROOFIT_LIBS)


# MyCB.o: MyCB.cc MyCB.h
# 	$(CXX) $(FLAGS) -c MyCB.cc -o MyCB.o $(ROOFIT_LIBS)

$(OBJ_DIR)/RooMcorMvisTsallis.o : $(SRC_DIR)/RooMcorMvisTsallis.cc $(INC_DIR)/RooMcorMvisTsallis.h
	$(CXX) $(FLAGS) -c $(SRC_DIR)/RooMcorMvisTsallis.cc -o $(OBJ_DIR)/RooMcorMvisTsallis.o $(ROOFIT_LIBS)

$(OBJ_DIR)/usefulFunctions.o: $(SRC_DIR)/usefulFunctions.cc $(INC_DIR)/usefulFunctions.h
	$(CXX) $(FLAGS) -c $(SRC_DIR)/usefulFunctions.cc -o $(OBJ_DIR)/usefulFunctions.o $(ROOFIT_LIBS)

$(OBJ_DIR)/fitter_utils.o: $(SRC_DIR)/fitter_utils.cc $(INC_DIR)/fitter_utils.h $(INC_DIR)/RooMcorMvisTsallis.h $(INC_DIR)/usefulFunctions.h
	$(CXX) $(FLAGS) -c $(SRC_DIR)/fitter_utils.cc -o $(OBJ_DIR)/fitter_utils.o $(ROOFIT_LIBS)

clean:
	rm -f $(OBJ_DIR)/*.o $(BIN_DIR)/*
