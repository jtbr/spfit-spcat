CC=gcc
CXX=g++
#CFLAGS=-g -Od -Wall -Wextra  # debug
#CFLAGS=-O3 -Wall  # optimized for distribution
CFLAGS=-O3 -Wall -Wextra -march=native -ffp-contract=off -Isrc  # optimized for speed on current device (-Ofast is faster but math results may differ very slightly)
CXXFLAGS=$(CFLAGS)
EXEQ=spfit spcat calmrg  # dpcat and dpfit are optional variations to spcat and spfit respectively
EXEA=${EXEQ} moiam stark termval sortn calbak reassign sortegy iambak iamcalc
EXAMPLES=fit_example cat_example
#next line for atlas blas
#BLASLIB=-lcblas -latlas
#OpenBLAS from libopenblas-dev apt package
#BLASLIB=-lopenblas
#if BLASLIB undefined, use supplied fallback BLAS routines
ifndef BLASLIB
	LBLAS=dblas.o
endif

vpath %.c   src/spfit src/spcat src/engine src/splib src/common src/legacy_apps
vpath %.cpp src/spfit src/spcat src/engine src/splib src/common src/legacy_apps examples
vpath %.h   src/spfit src/spcat src/engine src/splib src/common src/legacy_apps
vpath %.hpp src/spfit src/spcat src/engine src/splib src/common src/legacy_apps

default: ${EXEQ}
all: ${EXEA}
install:
	-mv ${EXEQ} /usr/local/bin
examples: ${EXAMPLES}
clean:
	rm -f ${EXEA} ${EXAMPLES} *.o *.a

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
fit_example: fit_example.o CalFit.o CalFit_helpers.o CalFitIO.o subfit.o ftran.o $(OBJFILES) $(LBLAS); g++ -o $@ $^ $(BLASLIB) -lm
cat_example: cat_example.o CalCat.o CalCat_helpers.o CalCatIO.o sortsub.o $(OBJFILES) $(LBLAS); g++ -o $@ $^ $(BLASLIB) -lm

splib.a: ulib.o cnjj.o catutil.o lsqfit.o $(LBLAS)
	ar r splib.a $^
	ranlib splib.a

fit_main.o: fit_main.cpp SpinvEngine.hpp DpiEngine.hpp CalFit.hpp CalFitIO.hpp lsqfit.h subfit.h CalError.hpp file_helpers.hpp SigintFlag.hpp
subfit.o: subfit.cpp calpgm_types.h blas_compat.h ulib.h subfit.h CalError.hpp file_helpers.hpp
lsqfit.o: lsqfit.c lsqfit.h cblas.h
cat_main.o: cat_main.cpp calpgm_types.h SpinvEngine.hpp DpiEngine.hpp CalCat.hpp CalCatIO.hpp OutputSink.hpp ulib.h sortsub.h CalError.hpp file_helpers.hpp SigintFlag.hpp
sortsub.o: sortsub.c calpgm_types.h sortsub.h
calmrg.o: calmrg.cpp calpgm_types.h ulib.h catutil.h CalError.hpp file_helpers.hpp SigintFlag.hpp
termval.o: termval.cpp calpgm_types.h ulib.h CalError.hpp file_helpers.hpp SigintFlag.hpp
readopt.o: readopt.c readopt.h calpgm_types.h ulib.h
stark.o: stark.cpp calpgm_types.h blas_compat.h ulib.h CalError.hpp file_helpers.hpp SigintFlag.hpp
iamcalc.o: iamcalc.cpp calpgm_types.h blas_compat.h ulib.h readopt.h CalError.hpp file_helpers.hpp SigintFlag.hpp
reassign.o: reassign.cpp calpgm_types.h ulib.h CalError.hpp file_helpers.hpp SigintFlag.hpp
moiam.o: moiam.cpp calpgm_types.h blas_compat.h ulib.h CalError.hpp file_helpers.hpp SigintFlag.hpp
ulib.o: ulib.c calpgm_types.h blas_compat.h ulib.h cnjj.h
cnjj.o: cnjj.c cnjj.h
file_helpers.o: file_helpers.cpp file_helpers.hpp CalError.hpp
SigintFlag.o: SigintFlag.cpp SigintFlag.hpp
Logger.o: Logger.cpp Logger.hpp
spinv_setup.o: spinv_setup.cpp calpgm_types.h blas_compat.h ulib.h cnjj.h spinit.h spinv_internal.h SpinvContext.hpp CalError.hpp
spinv_spin_symmetry.o: spinv_spin_symmetry.cpp calpgm_types.h blas_compat.h ulib.h cnjj.h spinit.h spinv_internal.h SpinvContext.hpp CalError.hpp
spinv_linalg_sort.o: spinv_linalg_sort.c calpgm_types.h blas_compat.h ulib.h cnjj.h spinit.h spinv_internal.h SpinvContext.hpp
spinv_hamiltonian.o: spinv_hamiltonian.cpp calpgm_types.h blas_compat.h ulib.h cnjj.h spinit.h spinv_internal.h SpinvContext.hpp CalError.hpp
spinv_utils.o: spinv_utils.cpp calpgm_types.h blas_compat.h ulib.h cnjj.h spinit.h spinv_internal.h SpinvContext.hpp CalError.hpp
spinit.o: spinit.cpp calpgm_types.h blas_compat.h ulib.h cnjj.h spinit.h spinv_internal.h SpinvContext.hpp CalError.hpp
SpinvEngine.o: SpinvEngine.cpp SpinvEngine.hpp SpinvContext.hpp spinv_internal.h
DpiEngine.o: DpiEngine.cpp DpiEngine.hpp DpiContext.hpp dpi.h
dpi.o: dpi.cpp dpi.h calpgm_types.h blas_compat.h ulib.h cnjj.h DpiContext.hpp CalError.hpp
ftran.o: ftran.c
sortn.o: sortn.cpp CalError.hpp file_helpers.hpp SigintFlag.hpp
dblas.o: dblas.c
calbak.o: calbak.cpp calpgm_types.h ulib.h catutil.h CalError.hpp file_helpers.hpp SigintFlag.hpp
sortegy.o: sortegy.cpp calpgm_types.h ulib.h CalError.hpp file_helpers.hpp SigintFlag.hpp
iambak.o: iambak.cpp calpgm_types.h ulib.h readopt.h CalError.hpp file_helpers.hpp SigintFlag.hpp
CalFit.o: CalFit.cpp CalFit.hpp CalFitIO.hpp CalculationEngine.hpp calpgm_types.h lsqfit.h ulib.h subfit.h CalError.hpp SigintFlag.hpp SpinvEngine.hpp DpiEngine.hpp
CalFit_helpers.o: CalFit_helpers.cpp CalFit.hpp CalFitIO.hpp CalculationEngine.hpp calpgm_types.h lsqfit.h ulib.h subfit.h CalError.hpp SpinvEngine.hpp DpiEngine.hpp
CalFitIO.o: CalFitIO.cpp CalFitIO.hpp CalFit.hpp CalculationEngine.hpp calpgm_types.h ulib.h subfit.h CalError.hpp
CalCat.o: CalCat.cpp CalCat.hpp OutputSink.hpp CalculationEngine.hpp calpgm_types.h blas_compat.h ulib.h catutil.h CalError.hpp SigintFlag.hpp
CalCat_helpers.o: CalCat_helpers.cpp CalCat.hpp OutputSink.hpp CalculationEngine.hpp calpgm_types.h blas_compat.h ulib.h CalError.hpp
CalCatIO.o: CalCatIO.cpp CalCatIO.hpp CalCat.hpp OutputSink.hpp CalculationEngine.hpp calpgm_types.h ulib.h CalError.hpp file_helpers.hpp
fit_example.o: fit_example.cpp SpinvEngine.hpp CalFit.hpp CalFitIO.hpp CalError.hpp Logger.hpp
cat_example.o: cat_example.cpp SpinvEngine.hpp CalCat.hpp CalCatIO.hpp OutputSink.hpp CalError.hpp Logger.hpp
