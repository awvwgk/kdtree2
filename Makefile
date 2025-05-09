BUILD_DIR := build
SRC_DIR := src
TEST_DIR := tests
LIB := $(BUILD_DIR)/libkdtree2.a

# Compiler detection
F90 ?= $(notdir $(shell command -v ifort 2>/dev/null || \
	command -v gfortran 2>/dev/null || echo "notfound"))

# Check if a valid compiler was found
ifeq ($(F90), notfound)
	$(error No Fortran compiler found.)
endif

# Compiler flags
ifeq ($(F90), ifort)
	FFLAGS := -warn all -O3 -qopenmp -fno-alias -module $(BUILD_DIR)
	ifeq ($(DEBUG), 1)
		FFLAGS += -g -check all -traceback
	endif
endif

ifeq ($(F90), gfortran)
	FFLAGS := -Wall -O3 -fopenmp -J $(BUILD_DIR)
	ifeq ($(DEBUG), 1)
		FFLAGS += -g -fcheck=all
	endif
endif

# All object files
OBJECTS := $(patsubst $(SRC_DIR)/%.f90, $(BUILD_DIR)/%.o, $(wildcard $(SRC_DIR)/*.f90))

# All tests
TESTS := $(patsubst $(TEST_DIR)/%.f90, $(TEST_DIR)/%, $(wildcard $(TEST_DIR)/*.f90))

# default rule
default: makedir all

.PHONY: all
all:	$(LIB) $(TESTS)

.PHONY: makedir
makedir:
	@mkdir -p $(BUILD_DIR)

.PHONY: clean
clean:
	$(RM) $(BUILD_DIR)/*.o $(BUILD_DIR)/*.mod $(LIB) $(TESTS)

$(LIB): $(OBJECTS)
	$(AR) rcs $@ $^

# How to get .o object files from .f90 source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(F90) -c -o $@ $< $(FFLAGS)

# How to build tests
$(TEST_DIR)/%: $(TEST_DIR)/%.f90
	$(F90) -o $@ $^ $(FFLAGS)

# Dependencies
$(TESTS): $(LIB)
