# Trans4D Makefile
GPP := g++
LIB_SUFFIX := dll
GPP_FLAGS := -std=c++11 -Iinclude -fPIC -Wall -Wextra
DEBUG_FLAGS := -g
BUILD_DIR := bin
OBJ_DIR := build
SRC_DIR := src

ifeq ($(OS),Windows_NT)
	LIB_SUFFIX = dll
else
	LIB_SUFFIX = so
endif

TRANS4D_LIB_FILES := $(SRC_DIR)/InitBd.cpp $(SRC_DIR)/InitEq.cpp $(SRC_DIR)/InitPs.cpp $(SRC_DIR)/InitVl.cpp
TRANS4D_LIB_FILES += $(SRC_DIR)/Trans4d.cpp $(SRC_DIR)/Trans4dCommon.cpp

all: trans4dlib trans4d

meson:
	meson $(OBJ_DIR)
	cd $(OBJ_DIR) && meson compile

trans4dlib:
	@mkdir -p $(BUILD_DIR)
	$(GPP) -shared $(GPP_FLAGS) $(DEBUG_FLAGS) $(TRANS4D_LIB_FILES) -o $(BUILD_DIR)/libtrans4d.$(LIB_SUFFIX)

trans4d:
	@mkdir -p $(BUILD_DIR)
	@cp src/*.txt $(BUILD_DIR)
	$(GPP) $(GPP_FLAGS) $(DEBUG_FLAGS) src/Trans4dExample.cpp -o $(BUILD_DIR)/trans4d.exe -Lbin -ltrans4d

clean:
	@rm -rf $(SRC_DIR)/*.o
	@rm -rf bin build
