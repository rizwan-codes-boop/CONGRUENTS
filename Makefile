SHELL := /bin/sh

# The published source uses GNU C nested functions, so Apple Clang is not a
# compatible compiler. A project-local conda-forge GNU compiler is vendored.
DEPENDENCY_ROOT ?= .
GCC_PACKAGE := $(firstword $(wildcard $(DEPENDENCY_ROOT)/.conda-pkgs/gcc_impl_osx-arm64-*/bin/arm64-apple-darwin*-gcc))
GCC_LIBDIR := $(firstword $(wildcard $(DEPENDENCY_ROOT)/.conda-pkgs/libgcc-devel_osx-arm64-*/lib/gcc/arm64-apple-darwin*/*))
CC := $(if $(GCC_PACKAGE),$(GCC_PACKAGE),gcc)

CPPFLAGS := -I. -ICR_spectra -I$(DEPENDENCY_ROOT)/.deps/include -I$(DEPENDENCY_ROOT)/vendor/cubature -DCONGRUENTS_USE_SYSTEM_OPENMP
CFLAGS := -std=gnu11 -O3 -fopenmp -Wall -Wextra \
          -Wno-unused-function -Wno-unused-parameter
LDLIBS := $(DEPENDENCY_ROOT)/.deps/lib/libgsl.a $(DEPENDENCY_ROOT)/.deps/lib/libgslcblas.a \
          $(DEPENDENCY_ROOT)/.deps/lib/libcubature.a -L$(GCC_LIBDIR) -L$(DEPENDENCY_ROOT)/.deps/lib \
          -lemutls_w -lm -Wl,-rpath,$(abspath $(DEPENDENCY_ROOT))/.deps/lib

BIN_DIR := bin
PROGRAMS := $(BIN_DIR)/create_interp_objects $(BIN_DIR)/spectra
MODEL_HEADERS := $(wildcard *.h CR_spectra/*.h)

.PHONY: all check run-precompute run clean compiler-check

all: compiler-check $(PROGRAMS)

compiler-check:
	@test -x "$(CC)" || { echo "GNU GCC is missing. See README.md."; exit 1; }
	@test -f $(DEPENDENCY_ROOT)/.deps/lib/libgsl.a -a -f $(DEPENDENCY_ROOT)/.deps/lib/libcubature.a || { echo "Local GSL/cubature libraries are missing. See README.md."; exit 1; }

$(BIN_DIR):
	mkdir -p $@

$(BIN_DIR)/create_interp_objects: create_interp_objects.c $(MODEL_HEADERS) | $(BIN_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LDLIBS) -o $@

$(BIN_DIR)/spectra: spectra.c $(MODEL_HEADERS) | $(BIN_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LDLIBS) -o $@

check: all
	@$(BIN_DIR)/create_interp_objects 2>&1 | grep -q '^Usage:'
	@$(BIN_DIR)/spectra 2>&1 | grep -q '^Usage:'
	@echo "Executable smoke checks passed."

run-precompute: all
	mkdir -p data
	OMP_NUM_THREADS=$${OMP_NUM_THREADS:-4} $(BIN_DIR)/create_interp_objects input/cat_nt.txt data

run: all
	mkdir -p data output output/tau_loss
	OMP_NUM_THREADS=$${OMP_NUM_THREADS:-4} $(BIN_DIR)/spectra input/cat_nt.txt data output

clean:
	rm -rf $(BIN_DIR)

# Python/C ABI 2 library. Existing production executables remain unchanged.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
LIB_SUFFIX := dylib
SHARED_FLAGS := -dynamiclib
else
LIB_SUFFIX := so
SHARED_FLAGS := -shared
endif
PYTHON ?= python3
SHARED_LIBRARY := build/libcongruents.$(LIB_SUFFIX)
comma := ,
NATIVE_LDFLAGS := $(if $(GCC_LIBDIR),-L$(GCC_LIBDIR)) $(if $(GCC_PACKAGE),-L$(DEPENDENCY_ROOT)/.deps/lib -Wl$(comma)-rpath$(comma)$(abspath $(DEPENDENCY_ROOT))/.deps/lib)
.PHONY: shared test-week1
shared: $(SHARED_LIBRARY)

build:
	mkdir -p $@

$(SHARED_LIBRARY): csrc/congruents.c csrc/preparation.c csrc/internal.h csrc/include/congruents.h $(MODEL_HEADERS) | build
	$(CC) -Icsrc/include $(CFLAGS) -fPIC -fvisibility=hidden $(SHARED_FLAGS) csrc/congruents.c csrc/preparation.c $(NATIVE_LDFLAGS) -lm -o $@

.PHONY: test-portable test
test-portable: shared
	PYTHONPATH=src $(PYTHON) -m unittest discover -s tests -p test_serial.py -v

test: shared build/direct_ionisation build/direct_preparation
	PYTHONPATH=src $(PYTHON) -m unittest discover -s tests -v

build/direct_ionisation: tests/direct_ionisation.c CR_spectra/ionisation.h physical_constants.h | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LDLIBS) -o $@

test-week1: shared build/direct_ionisation
	PYTHONPATH=src $(PYTHON) -m unittest discover -s tests -p test_week1.py -v

.PHONY: test-week2
build/direct_preparation: tests/direct_preparation.c $(MODEL_HEADERS) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LDLIBS) -o $@

test-week2: shared build/direct_preparation
	PYTHONPATH=src $(PYTHON) -m unittest discover -s tests -p test_week2.py -v
