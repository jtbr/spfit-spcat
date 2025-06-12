# Documentation for DPFIT and DPCAT

These programs which use the subroutine DPLC are used to calculate $^2\Pi$ energies and intensities with one nuclear spin using Hund’s case (a). These programs are no longer separate executables. Instead, use SPCAT and SPFIT with the `--dpi` command-line argument. The `dpcat` and `dpfit` shell scripts can be used as before, as they will do this for you (on Unix or Mac).

## Format of Quantum Numbers

Quantum numbers which are used in the files can be given in several formats:

Linear Pi States and [Q]

| Linear Pi States  | | | | | [Q]  |
|--------------|----------------|-------------|-----|-----|------|
| $J + 1/2$    | $\Omega + 1/2$ | $\lambda$   | $v$ | $F$ | [19] |
| $J + 1/2$    | $\Omega + 1/2$ | $\lambda$   | $v$ |  -  | [19] |
| $J + 1/2$    | $\Omega + 1/2$ | $\lambda$   | $F$ |  -  | [8]  |
| $J + 1/2$    | $\Omega + 1/2$ | $\lambda$   |  -  |  -  | [8]  |

The field QNFMT in the .cat file can be regarded as having 3 sub-fields: QFMT = $Q*100 + H*10 + NQN$, in which NQN is the number of quanta per state, H is a binary code the existence of half integer quanta for F, and Q is the number in square brackets in the table above. The least significant bit of H refers to the F quantum number and is 1 if F is half integer. The quantum number $\lambda$ has a sign for the inversion parity, for e levels the parity is $(-1)^{J+1/2}$, while for f levels the parity is opposite.

## FORMAT of the .lin file

**line 1-NLINE** [`12I3,freeform`]: QN,FREQ,ERR,WT

**QN** = 12-integer field of quantum numbers. Interpreted in a multiple I3 format as the quantum numbers for the line (upper quanta first, followed immediately by lower quanta). Unused fields can be used for annotation. The entire field is printed in file.fit

**FREQ** = frequency in MHz or wavenumbers

**ERR** = experimental error. Minus sign means that the frequency and error are in units of wavenumbers. FREQ and ERR will be converted internally to units of MHz.

**WT** = relative weight of line within a blend (normalized to unity by program)

**notes**: If an end-of-file is encountered before all the lines are read in, NLINE is set to the number read to that point. If successive lines have the same frequency, the lines will be treated as a blend and derivatives will be averaged using WT/ERR. Any lines with format errors will be ignored.
The freeform input begins in column 37 and extends to the end of the line. See the notes at the end of the next section for more on the freeform input.

## FORMAT of the .par file and .var files

**line 1**: title
**line 2** [`freeform`]: NPAR, NLINE, NITR, NXPAR, THRESH , ERRTST, FRAC, CAL
(only NPAR used by CALCAT)

**NPAR** = maximum number of parameters
**NLINE** = maximum number of lines
**sign NITR** : negative value in DPFIT enables line assignment diagnostics
**mag NITR** = maximum number of iterations
**NXPAR** = number of parameters to exclude from end of list when fitting special lines (see notes)
**THRESH** = initial Marquardt-Levenburg parameter
**ERRTST** = maximum [(obs-calc)/error]
**FRAC** = fractional importance of variance. Positive value means multiply parameter errors by FRAC. Negative value means multiply parameter errors by -FRAC * RMS * SQRT (NLINE / NDFREE), where NLINE is the number of blends and NDFREE = NLINE - (the number of free parameters).
**CAL** = scaling for infrared line frequencies

**line 3 option information**[`freeform`]: SPIND, NVIB

**SPIND** = degeneracy of nuclear spin
**NVIB** = number of vibrations

**Parameter lines** [`freeform`]: IDPAR, PAR, ERPAR / LABEL

**IDPAR** is a parameter identifier. If NVIB > 1, IDPAR = IV + 100 * IDPAR0, where IV is the vibrational or electronic quantum number. If NVIB = 1, IDPAR = IDPAR0. There are no matrix elements defined connecting the vibrational states. If the sign of IDPAR is negative, DPFIT constrains the ratio of this parameter to the previous parameter to a fixed value during the fit.

**PAR** is the parameter value
**ERPAR** is the parameter uncertainty
**LABEL** is a parameter label (up to 10 characters are used) that is delimited by /

PARAMETER identifiers (IDPAR0) are:

| IDPAR0 | Parameter                                         |
|--------|---------------------------------------------------|
| 1      | A                                                 |
| 2      | $A_J$                                             |
| 3      | $A_H$                                             |
| 4      | $B+q/2$                                           |
| 5      | D                                                 |
| 6      | H                                                 |
| 7      | p                                                 |
| 8      | q                                                 |
| 9      | $p_D$                                             |
| 10     | $q_D$                                             |
| 11     | $\chi_1 = 0.5 * [a - (b+c)/2]$                    |
| 12     | $\chi_2 = d/2$                                    |
| 13     | $\chi_3 = 1.5 * [a + (b+c)/2]$                    |
| 14     | $\chi_4$                                          |
| 15     | $\chi_{1D}$                                       |
| 16     | $\chi_{2D}$                                       |
| 17     | $\chi_{3D}$                                       |
| 18     | $\chi_{4D}$                                       |
| 19     | $\chi_{4Q}$                                       |
| 20     | $\zeta_1 = eQq_1 \text{ for } \Omega = 1/2$       |
| 21     | $\zeta_2 = eQq_1 \text{ for } \Omega = 3/2$       |
| 22     | $\zeta_3 = -0.5 * eQq_2$                          |
| 23     | $\zeta_{1D}$                                      |
| 24     | $\zeta_{2D}$                                      |
| 25     | $\zeta_{3D}$                                      |
| 26     | $\zeta_{3Q}$                                      |
| 27     | $\gamma$                                          |

**line (n+1)-end** [`8F10.6`]: ( ( V(i,j),j=1,i ),i=1,NPAR )
V = Choleski decomposition of the correlation matrix, optional for file.par (The index i is the index for the parameter.)
NOTE: Definitions for parameters 15-19 and 23-26 were changed on 9 July 2005. Subscript D is now half the anticommutator of the corresponding operator with $J(J+1) - S(S+1) = N(N +1) + 2 N \cdot S$. Subscript Q is now half the anticommutator of the corresponding (off-diagonal) operator with $\pm(J+1/2)$.

## FORMAT for the .int file

**line 1** [`freeform`]:
```
title
```
**line 2** [`freeform`]:
```
FLAGS,TAG,QROT,FBGN,FEND,STR0,STR1,FQLIM,TEMP
```

- **FLAGS** = IRFLG\*1000+OUTFLG\*100+STRFLG\*10+EGYFLG
  - **IRFLG** = 1 if constants are in wavenumbers
  - **IRFLG** = 0 if constants are in MHz
  - **OUTFLG** = 0 for short form file.out
  - **STRFLG** = 1 to enable file.str output
  - **EGYFLG** $\neq$ 0 to enable file.egy energy listing
  - **EGYFLG** = 2,4 to enable file.egy derivative listing
  - **EGYFLG** = 3,4 to enable file.egy eigenvector listing
  - **EGYFLG** = 4 to dump Hamiltonian with no diagonalization

- **TAG** = catalog species tag (integer)
- **QROT** = partition function for TEMP
- **FBGN** = beginning integer F quantum (round up)
- **FEND** = ending integer F quantum (round up)
- **STR0,STR1** = log strength cutoffs
- **FQLIM** = frequency limit in GHz
- **TEMP** = temperature for intensity calculation in degrees K (default is 300K)

**line 3-end** [`freeform`]:
```
IDIP,DIPOLE
```

- **IDIP** is coded in decimal digit form according to the format V2\*100+ V1
- **DIPOLE** = dipole value

## FORMAT of .cat **catalog output** file

[`F13.4,2F8.4,I2,F10.4,I3,I7,I4,12I2`]:
```
FREQ,ERR,LGINT,DR,ELO,GUP,TAG,QNFMT,QN
```

- **FREQ** = Frequency of the line
- **ERR** = Estimated or experimental error (999.9999 indicates error is larger)
- **LGINT** = Base 10 logarithm of the integrated intensity in units of nm$^2$ MHz
- **DR** = Degrees of freedom in the rotational partition function (0 for atoms, 2 for linear molecules, and 3 for nonlinear molecules)
- **ELO** = Lower state energy in wavenumbers
- **GUP** = Upper state degeneracy
- **TAG** = Species tag or molecular identifier. A negative value flags that the line frequency has been measured in the laboratory. The absolute value of TAG is then the species tag (as given in line 2 of file.int above) and ERR is the reported experimental error.
- **QNFMT** = Identifies the format of the quantum numbers given in the field QN.
- **QN(12)** = Quantum numbers coded according to QNFMT. Upper state quanta start in character 1. Lower state quanta start in character 14. Unused quanta are blank, quanta whose magnitude is larger than 99 or smaller than -9 are shown with alphabetic characters or \*\*. Quanta between -10 and -19 are shown as a0 through a9. Similarly, -20 is b0, etc., up to -259, which is shown as z9. Quanta between 100 and 109 are shown as A0 through A9. Similarly, 110 is B0, etc., up to 359, which is shown as Z9.

## Format of .str **strength output** file

 [`F15.4,E15.6,I5,1X,24A,I5`]:
 ```
FREQ,DIPOLE,QNFMT,QN,ITEM
```

- **FREQ** = Frequency of the line
- **DIPOLE** = Reduced matrix element of the transition dipole
- **QNFMT** = Identifies the format of the quantum numbers given in the field QN.
- **QN(12)** = Quantum numbers coded according to QNFMT. Upper state quanta start in character 1. Lower state quanta start in character 14. Unused quanta are blank, quanta whose magnitude is larger than 99 or smaller than -9 are shown with alphabetic characters or \*\*. Quanta between -10 and -19 are shown as a0 through a9. Similarly, -20 is b0, etc., up to -259, which is shown as z9. Quanta between 100 and 109 are shown as A0 through A9. Similarly, 110 is B0, etc., up to 359, which is shown as Z9.
- **ITEM** = identifies number of dipole

## Format of .egy **energy output** file

[`2I5,3F18.6,6I3`]:
```
IBLK,INDX,EGY,PMIX,ERR,QN
```

- **IBLK** = Internal Hamiltonian block number
- **INDX** = Internal index Hamiltonian block
- **EGY** = Energy in wavenumbers
- **ERR** = Expected error of the energy in wavenumbers
- **PMIX** = mixing coefficient
- **QN(6)** = Quantum numbers for the state
