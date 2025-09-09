#--------------------------------------------------------------------------
# Directories and file names

SRC_DIR = src
OBJ_DIR = build
MOD_DIR = src
DEP_FILE = dependencies.dep # Name of the dependencies file
EXE_FILE = ASTROBOULD # Name of the executable

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# Compilator configuration (Default)
ifndef INTEL
	INTEL = 0  
endif
ifndef AMD
	AMD = 0
endif
# Check if GCC variable defined
ifdef GCC
	INTEL = 0
	AMD = 0
endif


# Compilator
# Default to GNU Fortran Compiler
FC = gfortran

ifeq ($(INTEL),1)
# Use Intel Fortran Compiler if setvars is completed
FC = ifx
else ifeq ($(AMD),1)
# Use Amd Fortran Compiler if setvars is completed
FC = flang
endif


# Debug
ifdef DEBUG
MYLDFLAGS := -O0 -g
ifeq ($(INTEL),1)
MYFFLAGS := -traceback -fpe0 -fp-model=source -check all
else ifeq ($(AMD),1)
MYFFLAGS := -fbacktrace -ffpe-trap=zero,invalid,overflow,underflow -fpmodel=source -fcheck=bounds -fcheck=array-temps -fcheck=pointer
else
MYFFLAGS := -fcheck=all -fbacktrace -ffpe-trap=zero,invalid,overflow,underflow 
endif
else  # No debug
MYLDFLAGS := -O2
ifeq ($(INTEL),1)
MYFFLAGS := -fimf-domain-exclusion=15 -funsafe-math-optimizations -funroll-loops -finit-real=zero
else ifeq ($(AMD),1)
MYFFLAGS := -flang-experimental-exec -fapprox-func -fno-honor-infinities -fno-honor-nans
else
MYFFLAGS := -ffinite-math-only -funsafe-math-optimizations -funroll-loops -finit-real=zero
endif
endif

# Serial
ifdef PARALLEL
MYFFLAGS += -fopenmp
endif

# Compiler flags
ifeq ($(INTEL),1)
FFLAGS := -warn -march=x86-64-v3 -nogen-interfaces -Wshadow -std=f2008 $(MYFFLAGS)
else ifeq ($(AMD),1)
FFLAGS := -Wall -Wextra -march=x86-64-v3 -Wshadow -std=f2008 $(MYFFLAGS)
else
FFLAGS := -Wall -Wextra -march=native -Wshadow -std=f2008 $(MYFFLAGS)
endif
LDFLAGS = $(MYLDFLAGS)

#--------------------------------------------------------------------------

# Variables
## Source files
FIXED_SOURCES = $(wildcard $(SRC_DIR)/*.f)
FREE_SOURCES = $(wildcard $(SRC_DIR)/*.F90)
SRCS = $(FIXED_SOURCES) $(FREE_SOURCES)

## Objects
FIXED_OBJECTS = $(addprefix $(OBJ_DIR)/,$(notdir $(FIXED_SOURCES:.f=.o)))
FREE_OBJECTS = $(addprefix $(OBJ_DIR)/,$(notdir $(FREE_SOURCES:.F90=.o)))
OBJECTS = $(FIXED_OBJECTS) $(FREE_OBJECTS)

## Modules
MODULES = $(filter-out main.mod, $(notdir $(OBJECTS:.o=.mod)))

## Dependencies
MAKE_DEP_FILE = $(DEP_FILE)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# Rules

# Default rule
all: $(OBJ_DIR) main $(DEP_FILE)

# Pattern rules

# Parallel
parallel:
	$(MAKE) PARALLEL=1

# Debug
debug:
	$(MAKE) DEBUG=1

# Debug and parallel
debug_parallel:
	$(MAKE) DEBUG=1 PARALLEL=1


# Create build directory before compiling
$(OBJECTS): | $(OBJ_DIR)

# Create build directory
$(OBJ_DIR): 
	@mkdir -p $(OBJ_DIR)
	@{ \
	if [ $(INTEL) -eq 1 ]; then \
		echo "Compiling using: Intel Fortran Compiler"; \
	elif [ $(AMD) -eq 1 ]; then \
		echo "Compiling using: AMD Fortran Compiler"; \
	else \
		echo "Compiling using: GNU Fortran Compiler"; \
	fi;\
	}

# Compile source files .f to objects
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FFLAGS) -o $@ -c $< $(LDFLAGS)

# Compile source files .F90 to objects
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.F90
	$(FC) $(FFLAGS) -o $@ -c $< $(LDFLAGS)

# Link objects and create main
main: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXE_FILE) $^ $(LDFLAGS)
	@mkdir -p $(MOD_DIR)
	@mv -f $(MODULES) $(MOD_DIR) 2>/dev/null; true
	@touch $(EXE_FILE)
	@echo ""
	@echo "Compilation successful!"
	@if [ -n "$(PARALLEL)" ]; then echo " (Using OpenMP)"; fi
	@echo ""
	

# Create dependencies file with fortdepend
deps:
	@{ \
	if [ ! -f $(DEP_FILE) ]; then \
		echo "Creating dependencies file: $(DEP_FILE)"; \
		if command -v fortdepend >/dev/null 2>&1; then \
			fortdepend -w -o $(DEP_FILE) -f $(FREE_SOURCES) -b $(OBJ_DIR) -i omp_lib; \
		else \
			echo "Error: fortdepend is not installed. Please run 'make install' to install it."; \
			exit 1; \
		fi; \
	else \
		echo "Dependencies file already exists: $(DEP_FILE)"; \
		echo "Remove it and run 'make deps' to create a new one."; \
	fi; \
	}

install: # Install fortdepend
	@{ \
	if ! command -v fortdepend >/dev/null 2>&1; then \
		echo "Installing fortdepend"; \
		python -m pip install fortdepend; \
	else \
		echo "fortdepend already installed."; \
	fi; \
	}
	$(MAKE) deps

# Clean
clean:
	@echo "Cleaning up..."
	@echo "rm -rf *.o *.mod $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod $(EXE_FILE)"
	@rm -rf *.o *.mod $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod $(EXE_FILE)
	@{ \
	if [ -d $(OBJ_DIR) ] && [ -z "$$(ls -A $(OBJ_DIR))" ]; then \
		rmdir $(OBJ_DIR); \
		echo "rmdir $(OBJ_DIR)"; \
	fi; \
	}
	@{ \
	if [ -d $(MOD_DIR) ] && [ -z "$$(ls -A $(MOD_DIR))" ]; then \
		rmdir $(MOD_DIR); \
		echo "rmdir $(MOD_DIR)"; \
	fi; \
	}

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# Include (or not) and create (if necessary) main.dep
ifeq ($(filter-out main all debug,$(MAKECMDGOALS)),)
include $(DEP_FILE)
endif

${MAKE_DEP_FILE}:
	$(MAKE) deps

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# Phony rules
.PHONY: clean all deps install parallel debug debug_parallel

# Default rule
.DEFAULT_GOAL := all
