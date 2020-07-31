# #get OS name and set make variables
# UNAME := $(shell uname -s)
SHELL := powershell.exe
GPP := g++
DLLEXT := .dll
GPP_FLAGS := -std=c++11
DEBUG_FLAGS := -g

TRANS4D_LIB_FILES := initbd.cpp initeq.cpp initps.cpp initvl.cpp
TRANS4D_LIB_FILES += trans4d.cpp trans4dcommon.cpp

all: trans4dlib trans4d

trans4dlib:
	$(GPP) -shared $(GPP_FLAGS) $(DEBUG_FLAGS) $(TRANS4D_LIB_FILES) -o ../bin/libtrans4d$(DLLEXT)
	copy ./Data*.txt ../bin/

trans4d:
	$(GPP) $(GPP_FLAGS) $(DEBUG_FLAGS) trans4d_example.cpp -o ../bin/trans4d.exe -L../bin -ltrans4d

clean:
	echo $(DLLEXT)
	rm *.dll *.so *.o *.out *.gch