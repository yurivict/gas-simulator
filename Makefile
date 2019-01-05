##
## Copyright (C) 2019 by Yuri Victorovich. All rights reserved.
##

CXX=		clang++70
SRCS=		main.cpp
APP=		sim
IMG_LIB=	-lfreeimage
COPTS=		-Wall -Wno-unused-const-variable
CXXFLAGS=	-O3 -std=c++17 $(COPTS) -I/usr/local/include -L/usr/local/lib $(IMG_LIB)

$(APP): $(SRCS) Makefile
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(APP) $(SRCS)
