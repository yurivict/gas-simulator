
SRCS=		main.cpp
APP=		sim
IMG_LIB=	-lfreeimage
CXXFLAGS=	-O3 -std=c++17 -I/usr/local/include -L/usr/local/lib -Wall $(IMG_LIB)

$(APP): $(SRCS) Makefile
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(APP) $(SRCS)
