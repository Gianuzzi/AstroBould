FC = gfortran -g
MYFFLAGS = -ffinite-math-only -funsafe-math-optimizations -funroll-loops #-fcheck=all
FFLAGS = -Wall -Wextra -march=native $(MYFFLAGS)
MYLDFLAGS = -O2
LDFLAGS = $(MYLDFLAGS)

TARGETS = main
FIXED_SOURCES = $(wildcard *.f)
FREE_SOURCES = $(wildcard *.F90)
FIXED_OBJECTS = $(patsubst %.f,%.o,${FIXED_SOURCES})
FREE_OBJECTS = $(patsubst %.F90,%.o,${FREE_SOURCES})

MAKEDEP = deps
DEP_FILE_NAME = main.dep
DEP_FILE = $(DEP_FILE_NAME)

all: $(TARGETS) $(DEP_FILE_NAME)

${FIXED_OBJECTS} : %.o : %.f
	$(FC) $(FFLAGS) -o $@ -c $< $(LDFLAGS)

${FREE_OBJECTS} : %.o : %.F90
	$(FC) $(FFLAGS) -o $@ -c $< $(DEP_FILE_NAME) $(LDFLAGS)
	
%.o: %.F90
	$(FC) $(FFLAGS) -o $@ -c $< $(LDFLAGS)

main: ${FREE_OBJECTS} ${FIXED_OBJECTS}
	$(FC) $(FFLAGS) -o main $^ $(LDFLAGS)

deps:
	@echo "Creating dependencies file: $(DEP_FILE_NAME)"
	@{ \
	if command -v fortdepend >/dev/null 2>&1; then \
		fortdepend -w -o $(DEP_FILE_NAME) -f *.F90; \
	else \
		echo "Error: fortdepend is not installed. Please run 'make install' to install it."; \
		exit 1; \
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
	make deps

clean:
	@echo "rm -rf *.o *.mod main"
	@rm -rf *.o *.mod main

# Include (or not) and create (if necessary) main.dep
ifeq ($(filter-out main all,$(MAKECMDGOALS)),)
ifneq ("$(wildcard $(DEP_FILE_NAME))","")
include $(DEP_FILE_NAME)
else
include $(DEP_FILE)
endif
endif

## The creation of main.dep is made with the python package: fortdepend
## It can be installed via pip: 
### $ python -m pip install fortdepend
${DEP_FILE}:
	make deps

.PHONY: clean all deps install
