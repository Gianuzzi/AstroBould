#--------------------------------------------------------------------------
# Project configuration

SRC_DIR   = src
OBJ_DIR   = build
MOD_DIR   = src
DEP_FILE  = dependencies.dep
EXE_FILE  = ASTROBOULD

#--------------------------------------------------------------------------
# Compiler selection

# Defaults
FC       = gfortran
INTEL    ?= 0
AMD      ?= 0
OPT      ?= 2   # default optimization level if not specified

# Intel compiler
ifeq ($(INTEL),1)
  FC = ifx
endif

# AMD compiler (AOCC Flang)
ifeq ($(AMD),1)
  FC = flang
endif

#--------------------------------------------------------------------------
# Flags

# Base settings
STD        = -std=f2008
ARCH       = -march=native
WARN_GCC   = -Wall -Wextra -Wshadow
WARN_INTEL = -warn -nogen-interfaces -Wshadow
WARN_AMD   =    # (no good Fortran warning flags in AOCC yet)

#--------------------------------------------------------------------------
# Optimization / Debug

# Debug build
ifdef DEBUG
  MYCFLAGS = -O0 -g
  ifeq ($(INTEL),1)
    MYFFLAGS = -traceback -fpe0 -check all -fp-model=source
  else ifeq ($(AMD),1)
    MYFFLAGS = -fbacktrace -ffpe-trap=zero,invalid,overflow,underflow -fsanitize=address,undefined
  else
    MYFFLAGS = -fcheck=all -fbacktrace -ffpe-trap=zero,invalid,overflow,underflow -fsanitize=address,undefined
  endif
else
  # Release build
  MYCFLAGS = -O$(OPT)
  ifeq ($(INTEL),1)
    MYFFLAGS = -funsafe-math-optimizations -funroll-loops \
               -xHost -qopt-report -fimf-domain-exclusion=15
  else ifeq ($(AMD),1)
    MYFFLAGS = -ffast-math -funroll-loops -fvectorize
  else
    MYFFLAGS = -ffinite-math-only -funsafe-math-optimizations \
               -funroll-loops -ftree-vectorize -finit-real=zero
  endif
endif

#--------------------------------------------------------------------------
# Parallelization
ifdef PARALLEL
  MYFFLAGS += -fopenmp
endif

#--------------------------------------------------------------------------
# Final flags per compiler
ifeq ($(INTEL),1)
  FFLAGS  = $(WARN_INTEL) $(STD) $(ARCH) $(MYCFLAGS) $(MYFFLAGS)
else ifeq ($(AMD),1)
  FFLAGS  = $(WARN_AMD) $(STD) $(ARCH) $(MYCFLAGS) $(MYFFLAGS)
else
  FFLAGS  = $(WARN_GCC) $(STD) $(ARCH) $(MYCFLAGS) $(MYFFLAGS)
endif

LDFLAGS = $(MYCFLAGS) $(MYFFLAGS)

#--------------------------------------------------------------------------
# Sources, objects, modules

FIXED_SOURCES  = $(wildcard $(SRC_DIR)/*.f)
FREE_SOURCES   = $(wildcard $(SRC_DIR)/*.F90)
SRCS           = $(FIXED_SOURCES) $(FREE_SOURCES)

FIXED_OBJECTS  = $(addprefix $(OBJ_DIR)/,$(notdir $(FIXED_SOURCES:.f=.o)))
FREE_OBJECTS   = $(addprefix $(OBJ_DIR)/,$(notdir $(FREE_SOURCES:.F90=.o)))
OBJECTS        = $(FIXED_OBJECTS) $(FREE_OBJECTS)

MODULES        = $(filter-out main.mod,$(notdir $(OBJECTS:.o=.mod)))

#--------------------------------------------------------------------------
# Rules

all: $(OBJ_DIR) main $(DEP_FILE)

debug:
	$(MAKE) DEBUG=1

parallel:
	$(MAKE) PARALLEL=1

debug_parallel:
	$(MAKE) DEBUG=1 PARALLEL=1

intel:
	$(MAKE) INTEL=1

amd:
	$(MAKE) AMD=1

$(OBJ_DIR):
	@mkdir -p $@
	@echo "Compiler: $(FC)"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

main: $(OBJECTS)
	$(FC) -o $(EXE_FILE) $^ $(LDFLAGS)
	@mkdir -p $(MOD_DIR)
	@mv -f $(MODULES) $(MOD_DIR) 2>/dev/null || true
	@echo "Compilation successful!"
	@if [ -n "$(PARALLEL)" ]; then echo " (Using OpenMP)"; fi

deps:
	@if [ ! -f $(DEP_FILE) ]; then \
	  echo "Creating dependencies file: $(DEP_FILE)"; \
	  if command -v fortdepend >/dev/null 2>&1; then \
	    fortdepend -w -o $(DEP_FILE) -f $(FREE_SOURCES) -b $(OBJ_DIR) -i omp_lib; \
	  else \
	    echo "Error: fortdepend is not installed."; \
	    exit 1; \
	  fi; \
	else \
	  echo "Dependencies file already exists: $(DEP_FILE)"; \
	  echo "Remove it and run 'make deps' to recreate."; \
	fi

install:
	@if ! command -v fortdepend >/dev/null 2>&1; then \
	  echo "Installing fortdepend"; \
	  python -m pip install fortdepend; \
	else \
	  echo "fortdepend already installed."; \
	fi
	$(MAKE) deps

clean:
	@echo "Cleaning up..."
	@rm -rf *.o *.mod $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod $(EXE_FILE)
	@[ -d $(OBJ_DIR) ] && rmdir --ignore-fail-on-non-empty $(OBJ_DIR) || true
	@[ -d $(MOD_DIR) ] && rmdir --ignore-fail-on-non-empty $(MOD_DIR) || true

ifeq ($(filter-out main all debug parallel debug_parallel,$(MAKECMDGOALS)),)
-include $(DEP_FILE)
endif

.PHONY: clean all deps install parallel debug debug_parallel
.DEFAULT_GOAL := all
