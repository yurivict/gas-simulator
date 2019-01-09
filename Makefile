##
## Copyright (C) 2019 by Yuri Victorovich. All rights reserved.
##

#
# options
#

DBG_TRACK_PARTICLES=0 # enable the particle tracking feature
DBG_SAVE_IMAGES=0     # save images representing evolution of one particluar area

CXXFLAGS_SER=	-DUSE_PARALLELISM=0
CXXFLAGS_PAR=	-DNCPU=8 -DUSE_PARALLELISM=1
LDFLAGS_SER=
LDFLAGS_PAR=	-pthread

CXX=		clang++70
HEADERS=	ThreadPool.h
SRCS=		main.cpp
APP=		sim
APP_SER=	$(APP).ser
APP_PAR=	$(APP).par
IMG_LIB=	-lfreeimage
COPTS_INLINE=	#-mllvm -inline-threshold=100000
COPTS_BUILD=	-DDBG_TRACK_PARTICLES=$(DBG_TRACK_PARTICLES) -DDBG_SAVE_IMAGES=$(DBG_SAVE_IMAGES)
COPTS=		-Wall -Wno-unused-const-variable $(COPTS_INLINE) $(COPTS_BUILD)
COPTS+=		-Wno-return-std-move # for https://github.com/jarro2783/cxxopts/issues/160
CXXFLAGS+=	-O3 -std=c++17 $(COPTS) -I/usr/local/include -L/usr/local/lib $(IMG_LIB)

$(APP_PAR): $(SRCS) $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_PAR) $(LDFLAGS) $(LDFLAGS_PAR) -o $(APP_PAR) $(SRCS)

$(APP_SER): $(SRCS) $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_SER) $(LDFLAGS) $(LDFLAGS_SER) -o $(APP_SER) $(SRCS)

both: $(APP_PAR) $(APP_SER)
