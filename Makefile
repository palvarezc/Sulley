CXX = g++

BASEDIR = $(shell pwd)
LIBDIR = $(BASEDIR)/lib
EXEDIR = $(BASEDIR)/bin
BINDIR = $(BASEDIR)/bin
SRCDIR = $(BASEDIR)/src
OBJDIR = $(BASEDIR)/obj
INCDIR = $(BASEDIR)/include
OBJ_EXT=o



ROOT_FLAGS   = $(shell $(ROOTSYS)/bin/root-config --cflags) -D__ROOFIT_NOBANNER
ROOT_LIBS    = $(shell $(ROOTSYS)/bin/root-config --libs) -lTreePlayer -lTreeViewer -lMinuit -lXMLIO -lMLP -lRIO
#ROOT_LIBS = -L$(TMVALOC)/lib/ -lTMVA

FLAGS = $(ROOT_FLAGS)

FLAGS += -I$(INCDIR)
#USERLIBS = -L/home/hep/th1011/Documents/MVA/kfold/lib -lKFold
#USERLIBS += -L/home/hep/th1011/Documents/MVA/ukfold/lib -lUKFold
#USERLIBS += -L/home/hep/th1011/Documents/MVA/dendrology-master/lib -lDendrology

ROOT_LIBS += $(USERLIBS)
FLAGS += $(USERLIBS)

LCGBOOST = /usr/local
BOOST_INC_DIR  = -I$(LCGBOOST)/include/
BOOST_LIB_OPT  = -L$(LCGBOOST)/lib -lboost_filesystem

FLAGS += ${BOOST_LIB_OPT} ${BOOST_INC_DIR}

#FLAGS += -I/home/hep/th1011/Documents/MVA/kfold/include 
#FLAGS += -I/home/hep/th1011/Documents/MVA/ukfold/include 
#FLAGS += -I/home/hep/th1011/Documents/MVA/dendrology-master/include


ROOFIT_LIBS  = ${ROOT_LIBS} -lRooFitCore -lRooFit -lRooStats -lHistFactory -lTreePlayer -lTreeViewer
ROOFIT_LIBS += $(USERLIBS)


SRCS=$(filter-out $(wildcard $(SRCDIR)/_*), $(wildcard $(SRCDIR)/*.cpp))      #select all .cpp as source files, ignore the ones starting with "_"
EXES=$(filter-out $(wildcard $(SRCDIR)/_*), $(wildcard $(SRCDIR)/*.cc))       #select all .cc as "mains", ignore the ones starting with "_"
BINS=$(subst $(SRCDIR), $(EXEDIR),$(subst .cc,,$(EXES)))                      #create the executable names by removing the .cc, and puting them in the executable dir
BRSS=$(subst $(SRCDIR), $(EXEDIR),$(subst .cpp,,$(EXES)))
OBJS_HPP+=$(subst $(SRCDIR), $(OBJDIR),$(subst cpp,$(OBJ_EXT),$(SRCS)))       #create all the object file name by removing the .cpp and adding a .o
OBJS_EXE+=$(subst $(SRCDIR), $(OBJDIR),$(subst cc,$(OBJ_EXT),$(EXES)))
INCS+=$(subst $(SRCDIR), $(OBJDIR),$(subst cpp,$(OBJ_EXT),$(SRCS)))



all : $(SRCS) $(OBJS_HPP) $(OBJS_EXE) $(BINS)

$(SRCDIR)/%.cpp : $(INCDIR)/%.h

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(FLAGS) -c $^ -o $@  $(ROOFIT_LIBS) 

$(OBJDIR)/%.o : $(SRCDIR)/%.cc
	$(CXX) $(FLAGS) -c $^ -o $@ $(ROOFIT_LIBS)

$(BINDIR)/%: $(OBJS_HPP) $(OBJDIR)/%.o
	$(CXX) $(FLAGS) -o $@ -g $^ $(ROOFIT_LIBS)



.PHONY : clean

clean :
	rm -f $(OBJDIR)/*.o
