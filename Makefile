##
## Copyright (C) 2019 by Yuri Victorovich. All rights reserved.
##

#
# options
#

USE_PARALLELISM=1

ifeq ($(USE_PARALLELISM),1)
  #CXXFLAGS+=-fopenmp
  CXXFLAGS+=-DNCPU=8
  LDFLAGS+=-pthread
  EXT=.par
else
  EXT=.ser
endif

CXX=		clang++70
HEADERS=	ThreadPool.h
SRCS=		main.cpp
APP=		sim$(EXT)
IMG_LIB=	-lfreeimage
COPTS_INLINE=	#-mllvm -inline-threshold=100000
COPTS_BUILD=	-DUSE_PARALLELISM=$(USE_PARALLELISM)
COPTS=		-Wall -Wno-unused-const-variable $(COPTS_INLINE) $(COPTS_BUILD)
CXXFLAGS+=	-O3 -std=c++17 $(COPTS) -I/usr/local/include -L/usr/local/lib $(IMG_LIB)

$(APP): $(SRCS) $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(APP) $(SRCS)
