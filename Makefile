SHELL := /bin/sh

# The published source uses GNU C nested functions, so Apple Clang is not a
# compatible compiler. A project-local conda-forge GNU compiler is vendored.
GCC_PACKAGE := $(firstword $(wildcard .conda-pkgs/gcc_impl_osx-arm64-*/bin/arm64-apple-darwin*-gcc))
GCC_LIBDIR := $(firstword $(wildcard .conda-pkgs/libgcc-devel_osx-arm64-*/lib/gcc/arm64-apple-darwin*/*))
CC := $(GCC_PACKAGE)

CPPFLAGS := -I. -ICR_spectra -I.deps/include -Ivendor/cubature -DCONGRUENTS_USE_SYSTEM_OPENMP
CFLAGS := -std=gnu11 -O3 -fopenmp -Wall -Wextra \
          -Wno-unused-function -Wno-unused-parameter
LDLIBS := .deps/lib/libgsl.a .deps/lib/libgslcblas.a \
          .deps/lib/libcubature.a -L$(GCC_LIBDIR) -L.deps/lib \
          -lemutls_w -lm -Wl,-rpath,$(CURDIR)/.deps/lib

BIN_DIR := bin
PROGRAMS := $(BIN_DIR)/create_interp_objects $(BIN_DIR)/spectra
MODEL_HEADERS := $(wildcard *.h CR_spectra/*.h)

.PHONY: all check run-precompute run clean compiler-check

all: compiler-check $(PROGRAMS)

compiler-check:
	@test -x "$(CC)" || { echo "GNU GCC is missing. See README.md."; exit 1; }
	@test -f .deps/lib/libgsl.a -a -f .deps/lib/libcubature.a || { echo "Local GSL/cubature libraries are missing. See README.md."; exit 1; }

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
