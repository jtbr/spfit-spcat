CC=gcc
CXX=g++
#CFLAGS=-g -Od -Wall -Wextra  # debug
#CFLAGS=-O3 -Wall  # optimized for distribution
CFLAGS=-O3 -Wall -Wextra -march=native -ffp-contract=off -Isrc  # optimized for speed on current device (-Ofast is faster but math results may differ very slightly)
CXXFLAGS=$(CFLAGS)
EXEQ=spfit spcat calmrg  # dpcat and dpfit are optional variations to spcat and spfit respectively
EXEA=${EXEQ} moiam stark termval sortn calbak reassign sortegy iambak iamcalc
#next line for atlas blas
#BLASLIB=-lcblas -latlas
#OpenBLAS from libopenblas-dev apt package
#BLASLIB=-lopenblas
#if BLASLIB undefined, use supplied fallback BLAS routines
ifndef BLASLIB
	LBLAS=dblas.o
endif

vpath %.c   src/spfit src/spcat src/engine src/splib src/common src/legacy_apps
vpath %.cpp src/spfit src/spcat src/engine src/splib src/common src/legacy_apps
vpath %.h   src/spfit src/spcat src/engine src/splib src/common src/legacy_apps
vpath %.hpp src/spfit src/spcat src/engine src/splib src/common src/legacy_apps

default: ${EXEQ}
all: ${EXEA}
install:
	-mv ${EXEQ} /usr/local/bin
clean:
	rm -f ${EXEA} *.o *.a

SPLIB_OBJFILES=ulib.o cnjj.o catutil.o lsqfit.o
SPINV_OBJFILES=spinv_setup.o spinv_spin_symmetry.o spinv_linalg_sort.o spinv_hamiltonian.o spinv_utils.o
CPP_HELPERS=file_helpers.o SigintFlag.o Logger.o
OBJFILES=dpi.o spinit.o $(SPLIB_OBJFILES) $(SPINV_OBJFILES) SpinvEngine.o DpiEngine.o $(CPP_HELPERS)
spfit: fit_main.o CalFit.o CalFit_helpers.o CalFitIO.o subfit.o $(OBJFILES) $(LBLAS); g++ -o $@ $^ $(BLASLIB) -lm
spcat: cat_main.o CalCat.o CalCat_helpers.o CalCatIO.o sortsub.o $(OBJFILES) $(LBLAS); g++ -o $@ $^ $(BLASLIB) -lm
calmrg: calmrg.o $(CPP_HELPERS) splib.a; g++ -o $@ $^ $(BLASLIB) -lm
calbak: calbak.o $(CPP_HELPERS) splib.a; g++ -o $@ $^ $(BLASLIB) -lm
termval: termval.o $(CPP_HELPERS) splib.a; g++ -o $@ $^ $(BLASLIB) -lm
stark: stark.o $(CPP_HELPERS) splib.a; g++ -o $@ $^ $(BLASLIB) -lm
moiam: moiam.o $(CPP_HELPERS) ftran.o splib.a; g++ -o $@ $^ $(BLASLIB) -lm
iamcalc: iamcalc.o $(CPP_HELPERS) readopt.o ftran.o splib.a; g++ -o $@ $^ $(BLASLIB) -lm
#cnvwn: cnvwn.o splib.a ; gcc -o $@ $^
sortn: sortn.o $(CPP_HELPERS) sortsub.o; g++ -o $@ $^ -lm
reassign: reassign.o $(CPP_HELPERS) splib.a; g++ -o $@ $^ $(BLASLIB) -lm
sortegy: sortegy.o $(CPP_HELPERS) splib.a; g++ -o $@ $^ $(BLASLIB) -lm
iambak: iambak.o $(CPP_HELPERS) splib.a readopt.o; g++ -o $@ $^ $(BLASLIB) -lm

splib.a: ulib.o cnjj.o catutil.o lsqfit.o $(LBLAS)
	ar r splib.a $^
	ranlib splib.a

fit_main.o: fit_main.cpp calpgm.h SpinvEngine.hpp DpiEngine.hpp CalFit.hpp CalFitIO.hpp lsqfit.h
subfit.o: subfit.cpp calpgm.h
lsqfit.o: lsqfit.c lsqfit.h
cat_main.o: cat_main.cpp calpgm.h SpinvEngine.hpp DpiEngine.hpp CalCat.hpp CalCatIO.hpp
sortsub.o: sortsub.c calpgm.h
calmrg.o: calmrg.cpp calpgm.h CalError.hpp file_helpers.hpp SigintFlag.hpp
termval.o: termval.cpp calpgm.h CalError.hpp file_helpers.hpp SigintFlag.hpp
readopt.o: readopt.c readopt.h
stark.o: stark.cpp calpgm.h CalError.hpp file_helpers.hpp SigintFlag.hpp
iamcalc.o: iamcalc.cpp calpgm.h CalError.hpp file_helpers.hpp SigintFlag.hpp
reassign.o: reassign.cpp calpgm.h CalError.hpp file_helpers.hpp SigintFlag.hpp
moiam.o: moiam.cpp calpgm.h CalError.hpp file_helpers.hpp SigintFlag.hpp
ulib.o: ulib.c calpgm.h
cnjj.o: cnjj.c cnjj.h
file_helpers.o: file_helpers.cpp file_helpers.hpp CalError.hpp
SigintFlag.o: SigintFlag.cpp SigintFlag.hpp
Logger.o: Logger.cpp Logger.hpp
spinv_setup.o: spinv_setup.cpp calpgm.h spinit.h spinv_internal.h SpinvContext.hpp
spinv_spin_symmetry.o: spinv_spin_symmetry.cpp calpgm.h spinit.h spinv_internal.h SpinvContext.hpp
spinv_linalg_sort.o: spinv_linalg_sort.c calpgm.h spinit.h spinv_internal.h SpinvContext.hpp
spinv_hamiltonian.o: spinv_hamiltonian.cpp calpgm.h spinit.h spinv_internal.h SpinvContext.hpp
spinv_utils.o: spinv_utils.cpp calpgm.h spinit.h spinv_internal.h SpinvContext.hpp
spinit.o: spinit.cpp calpgm.h spinit.h spinv_internal.h SpinvContext.hpp
SpinvEngine.o: SpinvEngine.cpp SpinvEngine.hpp SpinvContext.hpp spinv_internal.h
DpiEngine.o: DpiEngine.cpp DpiEngine.hpp DpiContext.hpp dpi.h
dpi.o: dpi.cpp dpi.h calpgm.h DpiContext.hpp
ftran.o: ftran.c
sortn.o: sortn.cpp CalError.hpp file_helpers.hpp SigintFlag.hpp
dblas.o: dblas.c
calbak.o: calbak.cpp CalError.hpp file_helpers.hpp SigintFlag.hpp
sortegy.o: sortegy.cpp calpgm.h CalError.hpp file_helpers.hpp SigintFlag.hpp
iambak.o: iambak.cpp calpgm.h CalError.hpp file_helpers.hpp SigintFlag.hpp
CalFit.o: CalFit.cpp CalFit.hpp CalculationEngine.hpp
CalFit_helpers.o: CalFit_helpers.cpp CalFit.hpp
CalFitIO.o: CalFitIO.cpp CalFitIO.hpp CalFit.hpp CalculationEngine.hpp
CalCat.o: CalCat.cpp CalCat.hpp CalculationEngine.hpp
CalCat_helpers.o: CalCat_helpers.cpp CalCat.hpp
CalCatIO.o: CalCatIO.cpp CalCatIO.hpp CalCat.hpp CalculationEngine.hpp
