# Linking pythia to fastjet ( Adapted form Vincia1024 makefile and Dr. Soyez's tutorial)
# test OS - Ubuntu 11.04 64bit;
# Mushfiq, University of Dhaka , date : 06.11.2011 .

SHELL = /bin/sh
TOPDIR=$(shell \pwd)

# Need this to get SHAREDSUFFIX (e.g. dylib or so)
#-include ../config.mk

# -include ../.Makefile.inc

# deal with pythia 
#PYTHIADIR = its your local path to pythia. In my case its in the same directory
#where this Makefile is placed
PYTHIADIR   = /home/tjwilk/Physics/software/pythia8176
PYTHIAFLAGS = -I$(PYTHIADIR)/include 
PYTHIALIBS  =  -L$(PYTHIADIR)/lib/archive -lpythia8 -llhapdfdummy


# Include Pythia and Pythia examples config files
-include $(PYTHIADIR)/config.mk
-include $(PYTHIADIR)/examples/config.mk

# deal with fastjet
# FASTJETDIR = your local path to fastjet. 

FASTJETDIR  = /home/tjwilk/Physics/software/fastjet-install
FASTJETINC :=`$(FASTJETDIR)/bin/fastjet-config --cxxflags --plugins`
FASTJETLIB :=`$(FASTJETDIR)/bin/fastjet-config --libs --plugins`
CXX=g++
CXXFLAGS= -std=c++11 -O3 -Wall -Woverloaded-virtual -g -Wunused-parameter

# install_script = $(SHELL) ../utils/install-sh
# check_script = ../utils/check.sh


NAME=GlobalSoftDrop
NAME2=Nsubjettiness
SRCS=/home/tjwilk/Physics/thaler_urop/GlobalSoftDrop/GlobalSoftDrop.cc
SRCS+=/home/tjwilk/Physics/software/fjcontrib/Nsubjettiness/Nsubjettiness.cc 
# EXAMPLES=
INSTALLED_HEADERS=/home/tjwilk/Physics/thaler_urop/GlobalSoftDrop/GlobalSoftDrop.hh
INSTALLED_HEADERS+=/home/tjwilk/Physics/software/fjcontrib/Nsubjettiness/Nsubjettiness.hh /home/tjwilk/Physics/software/fjcontrib/Nsubjettiness/Njettiness.hh

# LHAPDF if you want to replace Pythia's dummy default
#LHAPDFARGS =
#LHAPDFARGS = -L$(HOME)/local/lib -lLHAPDF -I$(HOME)/local/include

# main61.cc in the example folder uses fastjet. Put this in the same folder
# where the make file is placed. 

# A few variables used in this Makefile:
# EX           := test_event_display
# EXE          := $(addsuffix .exe,$(EX))
#PYTHIA8      ?= $(PWD)/..
PYTHIA8	     ?= /home/tjwilk/Physics/software/pythia8176
STATICLIB    := $(PYTHIA8)/lib/archive/libpythia8.a
SHAREDLIB    := $(PYTHIA8)/lib/libpythia8.$(SHAREDSUFFIX)
DICTCXXFLAGS := -I$(PYTHIA8)/include
ROOTCXXFLAGS := $(DICTCXXFLAGS) $(shell root-config --cflags)
CXXFLAGS     += -I/home/tjwilk/Physics/software/fastjet-install/include -L$(PYTHIA8)/lib -lpythia8 -llhapdfdummy -L/home/tjwilk/Physics/software/fastjet-install/lib -lfastjet
# Libraries to include if GZIP support is enabled
ifeq (x$(ENABLEGZIP),xyes)
LIBGZIP=-L$(BOOSTLIBLOCATION) -lboost_iostreams -L$(ZLIBLOCATION) -lz
endif

# LDFLAGS1 for static library, LDFLAGS2 for shared library
LDFLAGS1 := $(shell root-config --ldflags --glibs) \
  -L$(PYTHIA8)/lib/archive -lpythia8 -llhapdfdummy $(LIBGZIP)
LDFLAGS2 := $(shell root-config --ldflags --glibs) \
  -L$(PYTHIA8)/lib -lpythia8 -llhapdfdummy $(LIBGZIP)

OBJS = $(SRCS:.cc=.o)
# EXAMPLES_SRCS = $(EXAMPLES:=.cc)

# install_HEADER  = $(install_script) -c -m 644
# install_LIB     = $(install_script) -c -m 644
# install_DIR     = $(install_script) -d
# install_DATA    = $(install_script) -c -m 644
# install_PROGRAM = $(install_script) -c -s
# install_SCRIPT  = $(install_script) -c
# Default target; make examples (but not shared dictionary)
#all: $(EX)

.PHONY: clean distclean examples check install

all: lib$(NAME).a 

lib$(NAME).a: $(OBJS) 
	ar cru lib$(NAME).a $(OBJS) 
	ranlib lib$(NAME).a

# examples: $(EXAMPLES)

# $(EXAMPLES): % : %.o all
# 	$(CXX) -o $@ $< -L. -l$(NAME) $(LDFLAGS)
higgstest:
	g++ -o GSD_higgs_test GSD_higgs_test.cc -std=c++11 $(STATICLIB) $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETINC) $(FASTJETLIB) $(DICTCXXFLAGS) $(ROOTCXXFLAGS) $(LDFLAGS1) $(INSTALLED_HEADERS) -L/home/tjwilk/Physics/thaler_urop/GlobalSoftDrop -l$(NAME) -L/home/tjwilk/software/fjcontrib/Nsubjettiness -l$(NAME2)

# ttbarstudy: 
# 	g++ -o ttbarstudy_njettiness ttbarstudy_njettiness.cc -std=c++11 $(STATICLIB) $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETINC) $(FASTJETLIB) $(DICTCXXFLAGS) $(ROOTCXXFLAGS) $(LDFLAGS1) $(INSTALLED_HEADERS) -l$(NAME) -l$(NAME2)

# zprimestudy: 
# 	g++ -o zprimestudy_njettiness zprimestudy_njettiness.cc -std=c++11 $(STATICLIB) $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETINC) $(FASTJETLIB) $(DICTCXXFLAGS) $(ROOTCXXFLAGS) $(LDFLAGS1) $(INSTALLED_HEADERS) -l$(NAME) -l$(NAME2)

# higgs2bbstudy: 
# 	g++ -o higgs2bbstudy_njettiness higgs2bbstudy_njettiness.cc -std=c++11 $(STATICLIB) $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETINC) $(FASTJETLIB) $(DICTCXXFLAGS) $(ROOTCXXFLAGS) $(LDFLAGS1) $(INSTALLED_HEADERS) -l$(NAME) -l$(NAME2)

# Rule to build zjets_sim example. Needs static PYTHIA 8 library
#qg_norm_test: $(STATICLIB) qg_norm_test.cc
	#$(CXX) $(ROOTCXXFLAGS) $@.cc -o $@.exe $(LDFLAGS1)

# Rule to build full dictionary
dict: $(SHAREDLIB)
	rootcint -f pythiaDict.cc -c $(DICTCXXFLAGS) \
           -DPYTHIA8_COMPLETE_ROOT_DICTIONARY \
           pythiaROOT.h pythiaLinkdef.h
	$(CXX) -shared -fPIC -o pythiaDict.$(SHAREDSUFFIX) pythiaDict.cc \
         -DPYTHIA8_COMPLETE_ROOT_DICTIONARY \
         $(ROOTCXXFLAGS) $(LDFLAGS2)


# Error messages if PYTHIA libraries don't exist
$(STATICLIB):
	@echo "Error: PYTHIA 8 archive library must be built"
	@false
$(SHAREDLIB):
	@echo "Error: PYTHIA 8 shared library must be built"
	@false

# Clean up
# clean:
# 	rm -f $(EXE) hist.root pythiaDict.* \
#                treeDict.cc treeDict.h pytree.root

clean:
	rm -f *~ *.o *.a

GlobalSoftDrop.o: GlobalSoftDrop.hh
# Nsubjettiness.o: Nsubjettiness.hh