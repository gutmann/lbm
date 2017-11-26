RM=rm
BUILD=build/
# Note, this preface is primarily aimed at figureing out which compiler to use if different environment variables aren't set
ifndef COMPILER
	ifndef F90
		ifdef FC
			F90=${FC}
		else
			F90=gfortran	# NOTE: if nothing is defined we default to gfortran
		endif
$(info F90 undefined, changed to ${F90})
	endif
	ifeq (${F90},gfortran)
		COMPILER=gnu
	endif
	ifeq (${F90},ifort)
		COMPILER=intel
	endif
endif

ifeq (${COMPILER}, gnu)
	ifndef F90
		ifdef FC
			F90=${FC}
		else
			F90=gfortran
		endif
$(info F90 undefined, changed to ${F90})
	endif

	FFLAGS=-Ofast -fimplicit-none -fopenmp # -O3
	LFLAGS=-lgomp
	MODOUTPUT=-J $(BUILD)
	ifeq (${MODE}, debug)
		FFLAGS=-g -fcheck=all -fbacktrace -fimplicit-none
	endif
endif

ifeq (${COMPILER}, intel)
	ifndef F90
		ifdef FC
			F90=${FC}
		else
			F90=ifort
		endif
$(info F90 undefined, changed to ${F90})
	endif

	FFLAGS=-O3 -xHost -u -qopenmp
	LFLAGS=-qopenmp
	MODOUTPUT=-module $(BUILD)
	ifeq (${MODE}, debug)
		FFLAGS=-debug -debug-parameters all -traceback -g -u -check all -check noarg_temp_created -CB
	endif
endif

LDFLAGS=-L${NETCDF}/lib ${NCAR_LIBS_NETCDF} ${LFLAGS}
FCFLAGS=${FFLAGS} -c -I${NCAR_INC_NETCDF} ${MODOUTPUT}

all: lbm_basic

clean:
	${RM} lbm_basic $(BUILD)*.o $(BUILD)*.mod

lbm_basic:$(BUILD)lbm_basic.o $(BUILD)io_routines.o
	${F90} ${LDFLAGS} $^ -o $@

$(BUILD)lbm_basic.o: $(BUILD)io_routines.o

$(BUILD)%.o:%.f90
	${F90} ${FCFLAGS} $< -o $@
