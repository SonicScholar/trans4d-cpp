# #get OS name and set make variables
# UNAME := $(shell uname -s)

GPP := g++
DLLEXT := .dll
GPP_FLAGS := -std=c++11
DEBUG_FLAGS := -g

TRANS4D_LIB_FILES := initbd.cpp Trans4dCommon.cpp

all: trans4dlib trans4d

trans4dlib:
	$(GPP) -shared $(GPP_FLAGS) $(DEBUG_FLAGS) $(TRANS4D_LIB_FILES) -o ../bin/libtrans4d$(DLLEXT)

trans4d:
	$(GPP) $(GPP_FLAGS) $(DEBUG_FLAGS) trans4d.cpp -o ../bin/trans4d.exe -L../bin -ltrans4d

clean:
	echo $(DLLEXT)
	rm *.dll *.so *.o *.out *.gch