##
## Copyright (C) 2019 by Yuri Victorovich. All rights reserved.
##

#
# options
#

USE_PARALLELISM=1
DBG_TRACK_PARTICLES=0 # enable the particle tracking feature
DBG_SAVE_IMAGES=0     # save images representing evolution of one particluar area

ifeq ($(USE_PARALLELISM),1)
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
COPTS_BUILD=	-DUSE_PARALLELISM=$(USE_PARALLELISM) -DDBG_TRACK_PARTICLES=$(DBG_TRACK_PARTICLES) -DDBG_SAVE_IMAGES=$(DBG_SAVE_IMAGES)
COPTS=		-Wall -Wno-unused-const-variable $(COPTS_INLINE) $(COPTS_BUILD)
COPTS+=		-Wno-return-std-move # for https://github.com/jarro2783/cxxopts/issues/160
CXXFLAGS+=	-O3 -std=c++17 $(COPTS) -I/usr/local/include -L/usr/local/lib $(IMG_LIB)

$(APP): $(SRCS) $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(APP) $(SRCS)

both:
	@$(MAKE) USE_PARALLELISM=0
	@$(MAKE) USE_PARALLELISM=1
