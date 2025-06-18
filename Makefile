CC=gcc
#CFLAGS=-g -Od -Wall -Wextra  # debug
#CFLAGS=-O3 -Wall  # optimized for distribution
CFLAGS=-O3 -Wall -Wextra -march=native -I.  # optimized for speed on current device (-Ofast is faster but math results may differ very slightly)
EXEQ=spfit spcat calmrg  # dpcat and dpfit are optional variations to spcat and spfit respectively
EXEA=${EXEQ} moiam stark termval sortn calbak reassign sortegy iambak # iamcalc is broken
#next line for atlas blas
#BLASLIB=-lcblas -latlas
#OpenBLAS from libopenblas-dev apt package
BLASLIB=-lopenblas
#next line for supplied fallback routines
#BLASLIB=
ifndef ($(BLASLIB))
	LBLAS=dblas.o
endif
default: ${EXEQ}
all: ${EXEA}
install:
	-mv ${EXEQ} /usr/local/bin
clean:
	rm -f ${EXEA} *.o *.a

SPLIB_OBJFILES=ulib.o cnjj.o slibgcc.o catutil.o lsqfit.o
SPINV_OBJFILES=spinv_setup.o spinv_spin_symmetry.o spinv_linalg_sort.o spinv_hamiltonian.o spinv_utils.o
OBJFILES=dpi.o spinit.o $(SPLIB_OBJFILES) $(SPINV_OBJFILES) SpinvEngine.o DpiEngine.o
spfit: fit_main.o CalFit.o CalFitIO.o subfit.o $(OBJFILES); g++ -o $@ $^ $(BLASLIB) -lm
spcat: cat_main.o sortsub.o $(OBJFILES); g++ -o $@ $^ $(BLASLIB) -lm
calmrg: calmrg.o splib.a; gcc -o $@ $^ $(BLASLIB) -lm
calbak: calbak.o splib.a; gcc -o $@ $^ $(BLASLIB) -lm
termval: termval.o splib.a; gcc -o $@ $^ $(BLASLIB) -lm
stark: stark.o splib.a ; gcc -o $@ $^ $(BLASLIB) -lm
moiam: moiam.o ftran.o splib.a; gcc -o $@ $^ $(BLASLIB) -lm
iamcalc: iamcalc.o ftran.o splib.a; gcc -o $@ $^ $(BLASLIB) -lm
#cnvwn: cnvwn.o splib.a ; gcc -o $@ $^ $(BLASLIB) -lm
sortn: sortn.o sortsub.o; gcc -o $@ $^
reassign: reassign.o splib.a ; gcc -o $@ $^ $(BLASLIB) -lm
sortegy: sortegy.o splib.a ; gcc -o $@ $^ $(BLASLIB) -lm
iambak: iambak.o splib.a readopt.o ; gcc -o $@ $^ $(BLASLIB) -lm

splib.a: ulib.o cnjj.o slibgcc.o catutil.o lsqfit.o $(LBLAS)
	ar r splib.a $^
	ranlib splib.a

# spinv.a: spinv_setup.o spinv_spin_symmetry.o spinv_linalg_sort.o spinv_hamiltonian.o spinv_utils.o
# 	ar r libspinv.a $^
# 	ranlib libspinv.a

fit_main.o:fit_main.cpp calpgm.h SpinvEngine.hpp DpiEngine.hpp CalFit.hpp CalFitIO.hpp lsqfit.h
subfit.o:subfit.c calpgm.h
lsqfit.o:lsqfit.c lsqfit.h
cat_main.o:cat_main.cpp calpgm.h DpiEngine.hpp
sortsub.o: sortsub.c calpgm.h
calmrg.o:calmrg.c calpgm.h
termval.o:termval.c calpgm.h
readopt.o:readopt.c readopt.h
stark.o:stark.c calpgm.h
iamcalc.o:iamcalc.c calpgm.h
reassign.o:reassign.c calpgm.h
moiam.o:moiam.c calpgm.h
ulib.o:ulib.c calpgm.h
cnjj.o:cnjj.c cnjj.h
slibgcc.o:slibgcc.c calpgm.h
spinv_setup.o:spinv_setup.c calpgm.h spinit.h spinv_internal.h SpinvContext.hpp
spinv_spin_symmetry.o:spinv_spin_symmetry.c calpgm.h spinit.h spinv_internal.h SpinvContext.hpp
spinv_linalg_sort.o:spinv_linalg_sort.c calpgm.h spinit.h spinv_internal.h SpinvContext.hpp
spinv_hamiltonian.o:spinv_hamiltonian.c calpgm.h spinit.h spinv_internal.h SpinvContext.hpp
spinv_utils.o:spinv_utils.c calpgm.h spinit.h spinv_internal.h SpinvContext.hpp
spinit.o:spinit.c calpgm.h spinit.h
SpinvEngine.o: SpinvEngine.cpp SpinvEngine.hpp SpinvContext.hpp spinv_internal.h
DpiEngine.o: DpiEngine.cpp DpiEngine.hpp DpiContext.hpp dpi.h
dpi.o:dpi.c dpi.h calpgm.h DpiContext.hpp
ftran.o:ftran.c
sortn.o: sortn.c
sortsub.o: sortsub.c
dblas.o: dblas.c
util.o:util.c
calbak.o:calbak.c
CalFit.o: CalFit.cpp CalFit.hpp CalculationEngine.hpp
CalFitIO.o: CalFitIO.cpp CalFitIO.hpp CalFit.hpp CalculationEngine.hpp
