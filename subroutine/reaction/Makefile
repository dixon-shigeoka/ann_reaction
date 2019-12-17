#===============================================================================
# Makefile
#===============================================================================

ISTDIR = ../../library
OBJDIR = ./lib
NAME   = mtsdriver.so
TARGET = $(ISTDIR)/$(NAME)
#FC = ifort
FC = gfortran
DEBUG =
DEBUGOMP =
FFLAGS =
DFLAGS =

ifeq ($(FC),ifort) #compile options for ifort
  #FFLAGS = -O2 -xHOST -module $(OBJDIR) -fpp
  FFLAGS = -O2  -module $(OBJDIR) -fpp/
  FFLAGS+= -ipo -inline-level=2 -inline-forceinline      # inline expansion
  FFLAGS+= -parallel -fopenmp                             # parallel calculation
  FFLAGS+= -shared -fPIC -lgomp
  #FFLAGS+= -module $(OBJDIR) -fpp
  #FFLAGS+= -simd # force vectorize at !DIR$ SIMD, may slow down in para cal.
  FFLAGS+= -mcmodel=large -shared-intel # use CPU memory of 2GB or more
  FFLAGS+= -traceback -g    # for debug
  #FFLAGS+= -traceback -g  -module $(OBJDIR) -fpp    # for debug
  #FFLAGS+= -warn all -check all -std -fpe0 -ftrapuv # for debug
  #FFLAGS+= -qopt-report=5 -qopt-report-phase=vec                              \
           -qopt-report-file=$(OBJDIR)/$(NAME).optrpt
  #FFLAGS+= -qopt-report=2 -qopt-report-phase=par                              \
           -qopt-report-file=$(OBJDIR)/$(NAME).optrpt
  DFLAGS = -gen-depformat=make -gen-depshow=nointr_mod
endif

ifeq ($(FC),gfortran) #compile options for gfortran
  FFLAGS = -O2 -ffree-form -fimplicit-none -J $(OBJDIR)
  FFLAGS+= -fopenmp                                      # parallel calculation
  FFLAGS+= -fbacktrace -g                                # for debug
  FFLAGS+= -Wall -fbounds-check -Wuninitialized                               \
           -ffpe-trap=invalid,zero,overflow               # for debug
  FFLAGS+= -shared -fPIC -lgomp
  DFLAGS = -cpp -MMD -MG
endif

ifeq ($(DEBUG),on) #make all DEBUG=on (command line)
  NAME   = mtsdriver.so
  TARGET = $(ISTDIR)/$(NAME)
#  FFLAGS= -traceback -g                                 # for debug
  FFLAGS+= -warn all -check all -std -fpe0 -ftrapuv
#  FFLAGS+= -fpp -module $(OBJDIR)
endif
ifeq ($(DEBUGOMP),on) #make all DEBUGOMP=on (command line)
  NAME   = mtsdriver.so
  TARGET = $(ISTDIR)/$(NAME)
  FFLAGS= -traceback -g                                 # for debug
  FFLAGS+= -parallel -qopenmp                            # parallel calculation
  FFLAGS+= -fpp -module $(OBJDIR)
endif

FOBJ+= mod_parameter.o         mod_globals.o                    \
       mod_function.o          imtss.o        \

.SUFFIXES : .o .f90

$(OBJDIR)/%.o : %.f90
	/bin/rm -f $(patsubst %.o,%.d,$@)
	$(FC) -c $(FFLAGS) $(DFLAGS) -MF $(patsubst %.o,%.d,$@) -o $@ $<


$(TARGET) : $(patsubst %,$(OBJDIR)/%,$(FOBJ))
	$(FC) $(FFLAGS) -o $@ $^

clean:
	/bin/rm -f $(OBJDIR)/*.o
	/bin/rm -f $(OBJDIR)/*.mod
	/bin/rm -f $(OBJDIR)/*__genmod.f90
	/bin/rm -f $(OBJDIR)/*.d
	/bin/rm -f $(TARGET)

commit: ##Commit&push branch to remote repository
	@bash .script/upload.sh

update: ##Update branch
	@bash .script/update.sh

mupdate: ## Update Master branch
	@bash .script/master_update.sh

# Generation code     >>>>
debugcode: ## generate mc2dgo_debug, mc2dgo_debugomp
	@bash .script/make_debugcode.sh

# Run validation >>>>

cmp: ## run debugcode in com/Validations
	@bash .script/runcmp.sh

manucmp: ## run debugcode in your select com/Validations
	@bash .script/runcmp.sh read


# for full auto validation >>>>

autocmp: debugcode cmp ## generate & run debugcode in com/Validations

all: clean $(TARGET)

.PHONY: clean all commit update mupdate
