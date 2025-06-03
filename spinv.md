# DOCUMENTATION for SPFIT and SPCAT

These programs use subroutines in SPINV.C to calculate energies and intensities for asymmetric rotors and linear molecules with up to 999 vibrational states and up to 9 spins. No distinction is made between electronic states and vibrational states, or between electronic and nuclear spins. SPFIT is used for fitting transitions and term values, with no requirement that the transitions obey any particular selection rules. SPFIT takes input files with extensions .par and .lin, copies the .par file to a .bak file, creates new text output files with extensions .par, .fit, .var. The .par and .var files follow essentially the same format and contain fitting parameters and optionally correlation information. The .fit file contains the results of the fit. SPCAT is used for predicting line positions and strengths. It takes the .var file as input along with an .int file that specifies limits for the calculation and contains the transition dipoles. The main output files for SPCAT use extensions .out and .cat, which are for general information and for the catalog output format, respectively. The .cat file follows the format of the JPL catalog, but does not have experimental data flagged. Auxiliary output files with extensions .egy and .str can also be requested. The .egy file can contain energies, derivatives with respect to the parameters, eigenvalues, and the undiagonalized Hamiltonian. The .str file contains a list of all transition dipole moments. The file names for SPFIT and SPCAT can be specified as command line arguments in any order. The first file name is used as the base file name for any files not explicitly specified. If no command line arguments are specified the program will give a prompt for the file names.

Some of the details of the program are described in H. M. Pickett, "The Fitting and Prediction of Vibration-Rotation Spectra with Spin Interactions," *J. Molec. Spectroscopy*, **148**, 371-377 (1991). The $I_{tot}$ spin-coupling basis is described in H. M. Pickett, 'Spin eigenfunctions and operators for the $D_n$ groups.' *J. Molec. Spectroscopy*, **228**, 659-663 (2004). The Euler series is described in H. M. Pickett, 'Use of Euler series to fit spectra with application to water.' *J. Molec. Spectroscopy*, **233**, 174-179 (2005).

## CONTENTS

- [Format of Quantum Numbers](#format-of-quantum-numbers)
- [Format of the .lin File](#format-of-the-lin-file)
- [Format of the .par and .var Files](#format-of-the-par-file-and-var-files)
- [Format of the .int File](#format-of-the-int-file)
- [Format of the .cat File](#format-of-cat-catalog-output-file)
- [Format of the .str File](#format-of-str-strength-output-file)
- [Format of the .egy File](#format-of-egy-energy-output-file)
- [Operator and Parameter Transformations](#operator-and-parameter-transformations)
- [Special Considerations for Linear Molecules](#special-considerations-for-linear-molecules)
- [Special Considerations for l-doubled States](#special-considerations-for-l-doubled-states)
- [Special Considerations for Symmetric Tops](#special-considerations-for-symmetric-tops)
- [Asymmetric Rotor Example](#asymmetric-rotor-example)

## Format of Quantum Numbers

Quantum numbers which are used in the files can be given in several formats:

| Linear Sigma States: |     |     |      |      |   | [Q]  |
| :------------------- | :-- | :-- | :--- | :--- | :-: | :--- |
| N                    | v   | J   | $F_1$  | $F_2$  | F | [12] |
| N                    | J   | $F_1$  | $F_2$  | $F_3$ | F | [01] |
| N                    | v   | J   | $F_1$  | $I_{tot}$ | F | [32] |
| N                    | v   | J   | ns  | $I_{tot}$ | F | [32]'|
| N                    | J   | $F_1$  | $F_2$  | $I_{tot}$ | F | [22] |
| N                    | J   | $F_1$  | $F_2$  | $I_{tot}$ | F | [21]'|
| N                    | J   | $F_1$  | ns  | $I_{tot}$ | F | [61] |
| N                    | v   | nn  | F   |     |    | [52] |
| N                    | nn  | F   |     |     |    | [41] |

| Symmetric Tops: |     |     |      |      |   | [Q]  |
| :-------------- | :-- | :-- | :--- | :--- | :-: | :--- |
| N               | K   | v   | J    | $F_1$  | F | [13] |
| N               | K   | J   | $F_1$  | $F_2$ | F | [02] |
| N               | K   | v   | J    | $I_{tot}$ | F | [33] |
| N               | K   | v   | ns   | $I_{tot}$ | F | [73] |
| N               | K   | J   | $F_1$ | $I_{tot}$ | F | [72] |
| N               | K   | J   | ns   | $I_{tot}$ | F | [62] |
| N               | K   | v   | nn   | F   |    | [53] |
| N               | K   |     | nn   | F   |    | [42] |

| **Asymmetric Tops:** |     |     |      |      |   |   [Q]   |
| :------------------- | :-- | :-- | :--- | :--- | :-: | :--- |
| N  | $K_a |  K_c$   | v   | J    | F | [14] |
| N  | $K_a |  K_c$   | J   | $F_1$  | F | [03] |
| N  | $K_a |  K_c$   | J   | $I_{tot}$ | F | [23] |
| N  | $K_a |  K_c$   | ns  | $I_{tot}$ | F | [63] |
| N  | $K_a |  K_c$   | v   | nn   | F | [54] |
| N  | $K_a |  K_c$   | nn  |  F   |   | [43] |

| **Hund's case $a$ Linear States:** |     |     |      |      |   |   [Q]   |
| :-------------- | :-- | :-- | :--- | :--- | :-: | :--- |
| J                    | $\Omega$ | $\Lambda$ | v    | $F_1$  | F | [14] |
| $J+1/2$              | $\Omega+1/2$ | $\Lambda$ | v    | $F_1$  | F | [19] |
| J                    | $\Omega$ | $\Lambda$ | $F_1$  | $F_2$ | F | [03] |
| $J+1/2$              | $\Omega+1/2$ | $\Lambda$ | $F_1$  | $F_2$ | F | [08] |

The quantum number **nn** is an aggregate spin quantum number which is used when the number of quantum numbers would otherwise be greater than 6. The quantum number **ns** is used as an auxiliary quantum number for n-fold symmetry with n = 3..6. Half integer spins are rounded up to the next integer. The sign of K for symmetric top notation designates parity under rotation about the b axis. When the vibronic wave function is even with respect to reflection in the ac plane, then the sign of K will also indicate the parity with respect to inversion. The symmetric top notation can also be used equivalently for linear molecules with orbital or vibrational angular momentum (lambda not zero). If the number of vibrations is one, then v is not included. The sequence $J, F_1... F$ can be replaced with $F_1, F_2... F$ if no electronic spin is present. The number of spins requested determines the length of the quantum number list. The factoring of the Hamiltonian is determined by the parameter set.

The field QNFMT in the .cat file can be regarded as having 3 sub-fields: QFMT = $Q*100 + H*10 + NQN$, in which NQN is the number of quanta per state, H is a binary code to indicate the existence of half integer quanta for the last three quantum numbers, and Q is the number in square brackets in the table above. (The least significant bit of H refers to the F quantum number and is 1 if F is half integer.)

The table entries for case a linear molecules are listed for completeness although SPCAT does not use these quanta. The entries with $\Omega + 1/2$ are for molecules with even electron spin multiplicity.

In most cases, the sequential spin coupling scheme is used, but an alternative $I_{tot}$ coupling that can be selected when there are n-equivalent spins.

| **Sequential Sum** |
| --- |
| $N + S = J$ |
| $J + I_1 = F_1$ |
| ... |
| $F_{m-2} + I_{m-1} = F_{m-1}$ |
| $F_{m-1} + I_m = F$ |

| **$I_{tot}$ Sum** |
| --- |
| $N + S = J$ |
| $J + I_1 = F_1$ |
| ... |
| $F_{m-n-1} + I_{m-n} = F_{m-n}$ |
| $F_{m-n} + I_{tot} = F$ |

where there are m nuclear spins, and $I_{tot}$ is the vector sum of the n-equivalent spins.
SPFIT and SPCAT are configured so that multiple spin types and symmetries can be handled as separate non-interacting vibronic states. The restrictions are that NQN and MOD(QNFMT / 100, 5) are the same for all vibronic states. The programs will append extra copies of F (corresponding to I = 0) if the number of spins is less than the maximum number of spins.

## FORMAT of the .lin file

**line 1-NLINE** [`12I3,freeform`]:
```
QN,FREQ,ERR,WT
```

- **QN** = 12-integer field of quantum numbers. Interpreted in a multiple I3 format as the quantum numbers for the line (upper quanta first, followed immediately by lower quanta). Unused fields can be used for annotation. The entire field is printed in file.fit
- **FREQ** = frequency in MHz or wavenumbers
- **ERR** = experimental error. Minus sign means that the frequency and error are in units of wavenumbers. FREQ and ERR will be converted internally to units of MHz.
- **WT** = relative weight of line within a blend (normalized to unity by program)

### Notes on line interpretation

If an end-of-file or empty line is encountered before all the lines are read in, NLINE is set to the number read to that point. If successive lines have the same frequency and experimental error, the lines will be treated as a blend and derivatives will be averaged using WT/ERR. Any lines with format errors will be ignored. If the blended lines are followed with a line with the same frequency and with an error that is at least 2 times larger than the error of the blend, then the error is interpreted as the uncertainty in the rms width of line frequencies within the blend, and the error-weighted square of the width is minimized along with the error-weighted square of the difference of experimental and calculated frequency. (The quantum numbers for this pseudo line are ignored.)

The freeform input begins in column 37 and extends to the end of the line. See the notes at the end of the next section for more on the freeform input.

## FORMAT of the .par file and .var files

**line 1** [`freeform`]:
```
title
```

**line 2** [`freeform`]:
```
NPAR, NLINE, NITR, NXPAR, THRESH , ERRTST, FRAC, CAL
```
(only NPAR is used by SPCAT)

- **NPAR** = maximum number of parameters
- **sign NLINE**: negative value allows up to 10 quanta per state
- **mag NLINE** = maximum number of lines
- **sign NITR** negative value in SPFIT enables line assignment diagnostics
- **mag NITR** = maximum number of iterations
- **NXPAR** = number of parameters to exclude from end of list when fitting special lines (see notes)
- **THRESH** = initial Marquardt-Levenburg parameter
- **ERRTST** = maximum [(obs-calc)/error]
- **FRAC** = fractional importance of variance. Positive value means multiply parameter errors by FRAC. Negative value means multiply parameter errors by -FRAC *RMS* SQRT (NLINE / NDFREE), where NLINE is the number of blends and NDFREE = NLINE - (the number of free parameters).
- **CAL** = scaling for infrared line frequencies

**line 3 option information** [`freeform`]:
```
CHR, SPIND, NVIB, KNMIN, KNMAX, IXX, IAX, WTPL, WTMN, VSYM, EWT, DIAG, XOPT
```
There can be one option line or multiple lines, as controlled by the option parameter VSYM. The first option line sets the default behavior of all the vibronic states, and successive lines (if present) modify the default behavior.

- **CHR** = character to modify parameter names file (must be in first column) sping.nam, default is 'g'. 'a' is used for Watson A set, 's' is used for Watson S set. Other character replaces the 'g' in the name 'sping'. Only used to label the .fit output file. (Ignored on all but first option line.) SPFIT looks for the .nam files in the current directory and then in the path given by the SPECNAME environment variable. (i.e. put something like SET SPECNAME=C:\SPECTRA\ in AUTOEXEC.BAT for Windows or setenv SPECNAME /spectra/ for unix). The trailing path delimiter is required.
- **sign SPIND** = If negative, use symmetric rotor quanta. If positive, use asymmetric rotor quanta (Sign ignored on all but first option line.)
- **mag SPIND** = degeneracy of spins, first spin degeneracy in units digit, second in tens digit, etc. (If last digit is zero, spin degeneracies occupy two decimal digits and the zero is ignored.) This field can specify up to 9 spins using up to 19 digits.
- **sign NVIB** = positive, prolate rotor (z = a, y = b, x = c) negative, oblate rotor (z = c, y = b, x = a) (Sign ignored on all but first option line.)
- **mag NVIB** = number of vibronic states on the first option line, identity of the vibronic state on all but the first option line. (max. value = 99 for 16-bit systems, max. value = 999 otherwise, vibrational quanta > 359 will be indicated by ** in the .cat file)
- **KNMIN, KNMAX** = minimum and maximum K values. If both = 0, then linear molecule is selected.
- **IXX** = binary flags for inclusion of interactions: bit 0 set means no $\Delta N \neq 0$ interactions, bit 1 means no delta J, bit 2 means no delta $F_1$, etc. [default = 0 includes all interactions] (Ignored on all but first option line.)
  - **mag IAX** = axis for statistical weight (1=a; 2=b; 3=c; 4= A, 2-fold top; 5=B, 2-fold top; 6= 3-fold top; 7=A, E, 4-fold top; 8=B, 4-fold top; 9=5-fold top; 10=A, E2, 6-fold top; 11=B, E1, 6-fold top). For mag IAX > 3, axis is b. (See Special Considerations for Symmetric Tops).
  - **sign IAX** = If negative, use $I_{tot}$ basis in which the last n spins are summed to give $I_{tot}$, which is then combined with the other spins to give F. (Sign is significant for all option lines.)
- **WTPL, WTMN** = statistical weights for even and odd state
- **VSYM** = If the value of VSYM on the first line is positive or zero, there is only one option line. The vibronic symmetry is coded as decimal digits (odd digit means reverse WTPL with WTMN) example: 10 = (v=0 even, v=1 odd) (Only works for the first 15 states) If VSYM is negative, signal that the next line is also an option line. If VSYM is positive or zero, signal that this is the last option line. In the multiple option line mode, the magnitude of VSYM is ignored.
- **sign EWT** If positive EWTFAC = 100. If negative EWTFAC = 1000. (Sign is ignored on all but first option line.)
- **EWT** = EWT0 + (EWT1 + EWT2) *EWTFAC = weight for states with E symmetry. (WTPL and WTMN apply to singly degenerate symmetries.) Ignore EWT0 if IAX < 6 or EWT0 = EWTFAC - 1. For the special case of n = 3 when WTPL = 0 or WTMN = 0, EWT0 should be doubled because only half of the states with MOD(K,3) $\neq$ 0 are calculated. See [Special Considerations for Symmetric Tops](#special-considerations-for-symmetric-tops). The weights WTPL, WTMN, and EWT0 can be divided by a common multiple if the rotational partition function is divided by the same factor.

    **Note**: If EWT1 = 0, then the state is non-degenerate with symmetry defined by IAX. Otherwise, the state is l-doubled and MUST be specified in adjacent pairs. The positive l state has EWT1 = l, and the negative l state has EWT1' = n' - |l|. For the pair of states, EWT1 + EWT1' = n'. For an n-fold top n' must be a multiple of n. When EWT1 = EWT1', the positive l is arbitrarily assigned to the lower v. States with MOD(EWT1, n) = 0 are $A_1$ and $A_2$ symmetry. States with MOD(EWT1,n) = n / 2 are $B_1$ and $B_2$ symmetry. Otherwise the states have E symmetry. E symmetry states will be designated with positive K for symmetric top quanta when n $\ge$ 3. For asymmetric top quanta with l > 0, $K_a + K_c = N+1$. For asymmetric top quanta with l < 0, $K_a + K_c = N$. If both WTPL and WTMN are not zero, there will be two states with the same nominal quantum number. (CALMRG will merge the degenerate transitions into a single line.) The value of EWT2 = 10 can be used to designate odd vibrational states.

- **DIAG** =
    -1 for no diagonalization
    0 for energy ordering within Wang sub-blocks
    1 for full projection assignment
    2 for energy ordering within Wang sub-blocks which follows order of diagonal elements of Hamiltonian
    3 for ordering by $\tau = K_a - K_c$ within vibration and spin sub-block set
    4 for ordering by $\sum < K, v|K_z^2|K, v > / \sum < K, v||K, v >$ within vibration and spin sub-block set
    5 for energy ordering within vibration and spin sub-block set which follows order of diagonal elements of Hamiltonian (Value ignored on all but first option line.)

- **XOPT** = PHASE + 10* NEWLZ + 20 *NOFC + 40* G12
    **PHASE** = 0 ... 8, force phase choice. 0 means no forcing, 8 means use standard phase with all even-order operators real and all odd-order operators imaginary. Other choices are a binary code for which even-order operators are imaginary.
    **NEWLZ** is 1 if Lz operator is explicitly defined in the parameter identifier (see [Special Considerations for $l$-doubled States](#special-considerations-for--doubled-states)).
    **NOFC** is 1 if the FF field of IDPAR specifies Kavg rather than Fourier coefficient.
    **G12** is one if the sine or cosine Fourier operator is to be multiplied by -1 if K is odd. (Ignored if NOFC is 1.)

For many cases only a single option line is needed. If different vibronic states have different spin multiplicity or different KMIN, KMAX additional lines are needed. Note that the sign of VSYM signals additional lines. The first option line sets up the defaults for all the vibrational states, and subsequent option lines specify deviations from the default. It is possible to mix Boson and Fermion states in the same calculation, e.g. fitting different isotopomers together, but the quantum number format (QNFMT) in SPCAT output will be correct only for the v=0 state.

**parameter lines** [`freeform`]:
```
IDPAR, PAR, ERPAR / LABEL
```

- **IDPAR** is a parameter identifier (see below). If the sign of IDPAR is negative, SPFIT constrains the ratio of this parameter to the previous parameter to a fixed value during the fit.
- **PAR** is the parameter value
- **ERPAR** is the parameter uncertainty
- **LABEL** is a parameter label (up to 10 characters are used) that is delimited by /

PARAMETER identifiers (IDPAR) are coded in decimal digit form in the order:

|              | EX | FF | I2 | I1 | NS | TYP | KSQ | NSQ | V2  | V1  |
| :----------- | :- | :- | :- | :- | :- | :-- | :-- | :-- | :-- | :-- |
| field digits | 1  | 2  | 1  | 1  | 1  | 2   | 1   | 1   | 1-3 | 1-3 |

- **EX** = Extended Euler-series Flag. (EX $\le$ 5) Operator is a term in an Euler series.
- **FF** = Fourier flag (used for internal rotation). If FF < 10, basic operator is multiplied by $\cos((FF) \cdot 2\pi K_{avg}\rho/3)$, where $\rho$ is coded by the absolute value of parameter ID=9100vv. See further discussion below. If $11 \le FF \le 20$, basic operator is multiplied by $\sin((FF-10) \cdot 2\pi K_{avg}\rho/3)$. If $21 \le FF \le 30$, basic operator is multiplied by $\cos((FF-10) \cdot 2\pi K_{avg}\rho/3)$, etc.
- **I2,I1** = spin identifiers [Usually I1 >= I2], I1=0 or I2=0 means N. If I1 = 0 and I2 > 0, perform commutator with $iN_z^2$.
- **NS** = power of $N \cdot S$ where S is the first spin. If NS > 4, subtract 5 and add $S_z N_z$ operator
- **TYP** = projection type
  - **0** = scalar
  - **1** = $N_a N_a$
  - **2** = $N_b N_b$
  - **3** = $N_c N_c$
  - **3+n** = $N_+^{2n} + N_-^{2n}$, n = 1 $\cdots$ 9 (L = 2n, $\Delta K = 2n$)
  - **11+n** = 'x' symmetry, n = 1 $\cdots$ 8 (L = 2n + 1, $\Delta K = 2n$)
  - **20+n** = off-diagonal 'a' symmetry, n = 0 $\cdots$ 19 (for prolate basis: L = n + 1, $\Delta K = 0, 2, 2, 4, 4, \ldots$)
  - **20** = $N_a$
  - **21** = $N_b N_c + N_c N_b$
  - **40+n** = off-diagonal 'b' symmetry, n = 0 $\cdots$ 19 (L = n + 1, $\Delta K = 1, 1, 3, 3, \ldots$)
  - **40** = $N_b$
  - **41** = $N_a N_c + N_c N_a$
  - **60+n** = off-diagonal 'c' symmetry, n = 0 $\cdots$ 19 (for prolate basis: L = n + 1, $\Delta K = 1, 1, 3, 3, \ldots$)
  - **60** = $N_c$
  - **61** = $N_a N_b + N_b N_a$
  - **80+n** = unique contribution for K' = K" = n $\cdot$ 10 + KSQ
  - **90+2n** = Euler series multiplying $N_+^{2n} + N_-^{2n}$
  - **91+2n** = constants for Euler and Fourier series
  - **91** = $\rho$ for Fourier series if KSQ = 0, NSQ1 = 0, and V1=V2. Only the absolute value is used for $\rho$, and the sign is used to designate a special symmetry (see below).
- **KSQ** = power of $N_z^2$
- **NSQ** = power of $N(N+1)$
- **V1, V2** = vibrational identifier. For NVIB < 9: V1 = V2 = 9 matches all V1 = V2. For 9 < NVIB < 99: V1 = V2 = 99 matches all V1 = V2. For 99 < NVIB < 360: V1 = V2 = 999 matches all V1 = V2. Usually, select V1 >= V2. V1 < V2 has special meaning for l-doubled states.

**Warning**: parameters requested really signify the operator that will multiply the parameter. Parameters not explicitly requested are presumed to be zero. For example, B for a linear molecule or symmetric top is 100, while for an asymmetric top it is 20000.

#### NOTES on parameter lines:

1. $N_{\pm} = N_x \pm iN_y$. If FF signifies a cosine series, then parameters with EVEN values of TYP > 20 and < 80 have an implicit i. If FF signifies a sine series, then parameters with ODD values of TYP > 20 and < 80 have an implicit i.
2. If EWT2 is the same for both of the vibrational $l$-doubled states **and** overall symmetry is $z$ or $y$, then the parameter is assumed to multiply $i L_z$ operator (see below). If EWT2 is different **and** overall symmetry is $x$ or *null*, then the parameter is also assumed to multiply $i L_z$ operator.
3. The sign of the $\rho$ parameter is used to designate a special symmetry for the Fourier series. If this sign is different for V1 and V2, then 0.5 is subtracted from NFF. For example, if NFF = 1, the basic operator is multiplied by $\cos(\pi K_{avg}\rho/3)$ instead of $\cos(2\pi K_{avg}\rho/3)$. If the magnitude of $\rho$ is not the same for the two states, replace $K_{avg}\rho$ with $(K_1\rho_1 + K_2\rho_2)/2$.
4. Prolate basis is $I^r$, and oblate basis is $III^l$. $\Delta K$ behavior for TYP = 20+n and 60+n are reversed for oblate basis.
5. $\rho_{vv}$ is specified by 9100vv. TYP=91+2n parameters, including $\rho$, are constants and are not fitted.
6. When NOFC = 0, the FF field is interpreted as a sample for a specific value of $K_{avg} = K_{avg}$ for $\Delta K$ even and $K_{avg} = K_{avg} + 1/2$ for $\Delta K$ odd. If FF < 10, basic operator is sampled at $K_{avg} = FF$. For l-doublets and a $\Delta l=0$ operator, the parameter is used for both $l$-doublets and the symmetry is the same as the corresponding cosine operator. If $11 \le FF \le 20$, basic operator is sampled at $K_{avg} = FF - 10$ and symmetry is the same as the corresponding sine operator. If $21 \le FF \le 30$, basic operator is sampled at $K_{avg} = FF -10$, etc.
7. For operators with I2 = 0 and I1 > 0, one value of N is replaced with the appropriate projection of $I_1$. For operators with I2 > 0 and I1 > 0, two values of N are replaced with appropriate projections of both I. For TYP=1,2,3, the operator is the expected Cartesian projection $I_{1z} I_{1z} / 3$. For example, 10000 is the $N_z N_z$ operator, 10010000 is the $N_z S_z$ operator, and 120010000 is the $S_z I_z - S \cdot I / 3$ operator.
8. For operators with I2 = I1 > 0, spin corrections appropriate for nuclear quadrupole coupling are applied: $SQRT((2I+1)/(2I-1))/4I$.
9. Whenever operators that do not commute are combined, the resulting operator is half the anti-commutator. The order of application is $N_z S_z$, followed by $N_z^2$, followed by $N^2$ and $N \cdot S$.
10. When coupling states where a given electronic (or nuclear) spin is different (e.g. spin-orbit coupling), the reduced matrix element for the spin operator, $< S||Op||S' >$, is assumed to be unity when $S > S'$ and $(-1)^{S'-S}$ when $S < S'$.
11. Euler denominator $a$ constants are defined by 9110vv', 9310vv',..., 9910vv'. Euler denominator $b$ constants are defined by 9101vv', 9301vv',..., 9901vv'. All operators with EX > 0 use the KSQ and NSQ fields to designate terms in an Euler series with denominator constants with TYP = $89 + 2 \cdot EX$.

**line (n+1)-end** [`8F10.6`]:
```
( ( V(i,j),j=1,i ) ,i=1,NPAR )
```
- **V** = Choleski decomposition of the correlation matrix, optional for file.par (The index i is the index for the parameter.)

### General Notes on .par and .var files

In the freeform input, the variables are all preset to reasonable default values. Any character not usually found in an E or F formatted number can separate the input numbers. A space or comma is recommended. Two successive commas indicate that the default value is to be used for that variable. At the end of the line or when a / character is encountered, all unspecified variables remain set to their default values. PAR defaults to zero. ERRPAR defaults to a very large number for SPFIT and to zero for SPCAT. If an end-of-file or error is encountered before the parameters are read in, NPAR is set to the number read to that point. If an end-of-file or error is encountered before V is completely read in, V is set to a unit matrix. The variables that are decoded to integer fields (such as IDPAR) cannot be larger than the number of significant digits in a double precision number (typically 15 digits).

Special lines to which NXPAR applies are lines in which the F quantum number is negative. In the quantum number assignment process in the program, the line is flagged and F is set to an appropriate value. When derivatives are accumulated,
the last NXPAR derivatives are ignored, and the energies are corrected by subtracting the first order contribution of these
parameters. If F < -1, the absolute value of F is used in the energy calculation. If F = -1, the F used is as close to the
previous spin quantum number as angular momentum addition rules allow. The value selected will be shown in the fit file
line listing in place of the -1.

If IDPAR is less than zero the magnitude is taken. In CALFIT, the parameter value will be constrained to be a constant
ratio of the preceding parameter value. In this way linear combinations of parameters can be fit as a unit.

Further discussion of the transformations between different sets of operators is given in [Operator and Parameter Trans-
formations](operator-and-parameter-transformations).


## FORMAT of the .int file

**line 1** [`freeform`]:
```
title
```

**line 2** [`freeform`]:
```
FLAGS,TAG,QROT,FBGN,FEND,STR0,STR1,FQLIM,TEMP,MAXV
```

- **FLAGS** = IRFLG\*1000+OUTFLG\*100+STRFLG\*10+EGYFLG
  - **IRFLG** = 1 if constants are in wavenumbers
  - **IRFLG** = 0 if constants are in MHz
  - **OUTFLG** = 0 for short form file.out
  - **STRFLG** = 1,2 to enable file.str output. SRFLG = 2 is used to label separate dipole contributions in the file.str output.
  - **EGYFLG** $\neq$ 0 to enable file.egy energy listing
  - **EGYFLG** = 2,4 to enable file.egy derivative listing
  - **EGYFLG** = 3,4 to enable file.egy eigenvector listing
  - **EGYFLG** = 5 to dump Hamiltonian with no diagonalization
- **TAG** = catalog species tag (integer)
- **QROT** = partition function for TEMP
- **FBGN** = beginning integer F quantum (round up)
- **FEND** = ending integer F quantum (round up)
- **STR0,STR1** = log strength cutoffs
- **FQLIM** = frequency limit in GHz
- **TEMP** = temperature for intensity calculation in degrees K (default is 300K)
- **MAXV** states with v>MAXV will not be included in the output files of SPCAT default is 999

**line 3-end** [`freeform`]:
```
IDIP,DIPOLE
```

- **IDIP** = dipole identifier (see below). A negative sign means that if STRFLG = 2 above then this dipole will be grouped with the previous dipole with IDIP positive.
- **DIPOLE** = dipole value

Dipole identifiers (IDIP) are coded in decimal digit form in the order:

|              | FC | TYP | I1 | V2  | V1  | SYM |
| :----------- | :- | :-- | :- | :-- | :-- | :-- |
| field digits | 1  | 1   | 1  | 1-3 | 1-3 | 1   |

- **FC** = second digit of Fourier order for TYP = 7,8; second digit of TYP otherwise
- **TYP** = dipole type
- **I1** = spin identifier [I1 = 0 means N or null]
- **V1,V2** = vibrational states [V1 usually >= V2]. Number of digits follows IDPAR in the .par file.
- **SYM** = symmetry [0 = magnetic, 1 = a, 2 = b, 3 = c]

| TYP | SYM  | absolute value of $\Delta K$ | Description                 |
| :-- | :--- | :----------- | :---------------------------------------------------------------------------------------------------------------------------------------- |
| 0   | 0    | 0            | magnetic dipole [N, S, I]                                                                                                                 |
| 0   | 1    | 0            | a dipole, $\phi_a$ if I1=0                                                                                                                |
| 0   | 2    | 1            | b dipole, $\phi_b$ if I1=0                                                                                                                |
| 0   | 3    | 1            | c dipole, $\phi_c$ if I1=0                                                                                                                |
| 1   | 0    | 0            | $(2\phi_z N_z - \phi_x N_x - \phi_y N_y)/2$ or $(2\phi_z I_z - \phi_x I_x - \phi_y I_y)/2$                                                    |
| 1   | 1    | 2            | $i (\{\phi_y, N_x\} + \{\phi_x, N_y\}) /2$ or $i (\{\phi_y, I_x\} + \{\phi_x, I_y\}) /2$                                                     |
| 1   | 2    | 1            | $i (\{\phi_z, N_y\} + \{\phi_y, N_z\}) /2$ or $i (\{\phi_z, I_y\} + \{\phi_y, I_z\}) /2$                                                     |
| 1   | 3    | 1            | $i (\{\phi_x, N_z\} + \{\phi_z, N_x\}) /2$ or $i (\{\phi_x, I_z\} + \{\phi_z, I_x\}) /2$                                                     |
| 2   | 0    | 2            | $\phi_y N_y - \phi_x N_x$ or $\phi_y I_y - \phi_x I_x$                                                                                      |
| 2   | 1-3  | TYP=1        | same as TYP = 1                                                                                                                           |
| 3   | any  | TYP=0        | $\{N_z^2, \text{TYP = 0 operator} \}/2$                                                                                                   |
| 4   | any  | TYP=0        | $\{N_z^2, \text{TYP = 0 operator} \}/2$                                                                                                   |
| 5   | 0    | 0            | $L_z \phi_z$                                                                                                                              |
| 5   | 1    | 0            | $[N_z^2, \phi_z]/2 = i (\{\phi_y, N_x\} - \{\phi_x, N_y\})/2$                                                                               |
| 5   | 2    | 1            | $[N_z^2, \phi_x]/2 = i (\{\phi_z, N_y\} - \{\phi_y, N_z\}) /2$                                                                               |
| 5   | 3    | 1            | $[N_z^2, \phi_y]/2 = i (\{\phi_x, N_z\} - \{\phi_z, N_x\}) /2$                                                                               |
| 6   | 0    | 2            | magnetic dipole [N, S, I] $\times (N_x^2 - N_y^2)$                                                                                        |
| 6   | 1    | 2            | $(\{\phi_z, N_x^2 - N_y^2\} + \{\phi_y, \{N_y, N_z\}\} - \{\phi_x, \{N_x, N_z\}\}) /\sqrt{6}$                                               |
| 6   | 2    | 3            | $(\{\phi_y, N_z^2 - N_x^2\} - \{\phi_z, \{N_x, N_y\}\}) /2$                                                                                 |
| 6   | 3    | 3            | $(\{\phi_x, N_y^2 - N_z^2\} + \{\phi_y, \{N_z, N_x\}\}) /2$                                                                                 |
| 7   | any  | TYP=0        | TYP=0 operator $\times \cos((FC*10 + I1) \cdot 2\pi K_{avg}\rho/3)$                                                                          |
| 8   | any  | TYP=0        | TYP=0 operator $\times i\sin((FC*10 + I1) \cdot 2\pi K_{avg}\rho/3)$                                                                         |
| 9   | any  | TYP=0        | $\{N_z^4, \text{TYP = 0 operator} \}/2$                                                                                                   |
| 10  | 1-3  | TYP=0        | $[N_z^2, [N_z^2, \text{TYP = 0 operator} ]]/4$                                                                                             |
| 11  | any  | TYP=2        | $[N_z^2, \text{TYP = 2 operator} ]/2$                                                                                                     |
| 12  | any  | TYP=2        | $i \{N_z, \text{TYP = 2 operator} \}/2$                                                                                                    |

**Notes**:

1. SYM = 1,2,3 maps to a,b,c regardless of whether the basis is oblate or prolate (as specified in the .par or .var file). Dipoles with SYM = 1,2,3 are assumed to be in units of Debye. Dipoles with SYM = 0 are assumed to be in units of a Bohr magneton. In the descriptions in the table above, [abc] = [zxy] for a prolate basis and abc = xyz for an oblate basis.
2. The column labeled $|\Delta K|$ gives the change in $K_a$ for a prolate basis or $K_c$ for an oblate basis.
3. Dipoles which are even order in direction cosine or N are assumed to be imaginary, except between l-doubled states with EWT1 $\neq$ 0.
4. The vector N for TYP = 0 (SYM = 0) is equivalent to $\phi_x N_x + \phi_y N_y + \phi_z N_z$
5. When the symmetry is 3-fold or lower (IAX $\le$ 6), dipoles between states with EWT1= (0,2), (2,0), and (2,2) are ignored, but the matrix elements are calculated using corresponding dipoles from states with EWT1=1.
6. For TYP = 7 or TYP = 8, I1 is used for the Fourier order and not the spin type. The constant $\rho$ is specified in the parameter set. The sign of the $\rho$ parameter is used to designate a special symmetry for the Fourier series. If this sign is different for V1 and V2, then 0.5 is subtracted from the Fourier order. For example, if IDIP = 72012, the basic b-dipole operator is multiplied by $\cos(3\pi K_{avg}\rho/3)$ instead of $\cos(4\pi K_{avg}\rho/3)$. If the magnitude of $\rho$ is not the same for the two states, replace $K_{avg}\rho$ with $(K_1\rho_1 + K_2\rho_2)/2$. TYP = 8 (with I1 > 0) dipoles are multiplied by i, and the symmetry of the states connected is 3 - SYM and the units follow the state symmetry (e.g. 81000 is in Debye). If the flag NOFC is set in the third line of the .var file, the dipole is sampled at fixed values of $K_{avg}$ as is done for IDPAR in the .par and.var files.
7. For ITYP = 7 or ITYP = 8, the Fourier order is I1 + 10 * FC. (FC $\le$ 9).
8. TYP = 1, 2, 5 are used for first-order Herman-Wallis corrections. TYP = 3, 4, 6, 10, 11, 12 are used for second-order Herman-Wallis corrections. Operators (TYP = 5,11) that involve a single commutator with $N^2$ have zero-valued matrix elements for Q-branch connections and have opposite signs for P and R branches. Operators (TYP = 10) that involve two commutators with $N^2$ also have zero-valued matrix elements for Q-branch connections but have the same signs for P and R branches. Further discussion of the transformations between different sets of dipole operators is given in Operator and Parameter Transformations.

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
- **DIPOLE** = Reduced matrix element of the transition dipole. If dipoles were magnetic, matrix element is multiplied by 0.009274 Debye / Bohr magneton.
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

## Operator and Parameter Transformations

The Hamiltonian or the transition dipole can be considered to be a linear dot product of parameters and operators. Two sets of parameters and operators are related by

$H = p^T Q = p'^T Q'$

where p is a parameter and Q is an operator. Let

$Q' = MQ$

where M is the transformation matrix. It is useful if M is square and not singular. Then

$H = p^T Q = p'^T Q' = p^T M^{-1} M Q$

and

$p' = (M^{-1})^T p \quad \text{or} \quad p = M^T p'$

Note that $(M^{-1})^T$ is only equal to M when M is orthogonal. In the more general case, the transformation of parameters and operators is not the same.
As a first example of these formal transformations, consider construction of a different operator, Q' from the standard SPCAT operator set, Q, using the constraints imposed by using a negative IDPAR in the .par file (or negative IDIP in the .int file). The constraint coefficients are a row in M. Related operators or additional composite operators are required to fill out M if we wish to determine p' in terms of p. This is also why changes in one element of Q' can have an impact on many other elements of p'.
Alternativly, Q' can represent a set of SPCAT operators, and Q can be a set of operators from the literature. As an example consider the second-order Hermann-Wallis operators for a b dipole. Eq. (2) can be written as

$Q' = \begin{pmatrix} D_0 \\ D_3 \\ D_4 \\ D_6 \\ D_{10} \\ D_{11} \end{pmatrix} = MQ = \begin{pmatrix} 1 & 0 & 0 & 0 & 0 & 0 \\ 0 & 1 & 1 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 & 0 & 0 \\ 0 & -1 & 1 & 0 & 2 & 0 \\ -1/2 & 0 & 1 & 1 & -1 & -1 \\ 0 & 0 & 1 & -1 & -1 & 1 \end{pmatrix} \begin{pmatrix} \phi_x \\ \phi_x N_z^2 \\ \{\phi_x, N_y^2\}/2 \\ \{\phi_x, N_x^2\}/2 \\ \{\phi_y, \{N_x,N_y\}\}/4 \\ \{\phi_z, \{N_x,N_z\}\}/4 \end{pmatrix}$

where the subscripts for D are the values of TYP in the .int file, and $\Phi$ is a standard cartesian set of operators. Then the parameters for the .int file, p', can be obtained from p using $(M^{-1})^T$. However, p, can be obtained from p' using only $M^T$ without inversion of M.
In this example,

$(M^{-1})^T = \begin{pmatrix} 1 & -1/8 & 1/8 & 0 & -1/8 & -1/4 \\ 0 & 3/4 & 1/4 & 0 & 1/4 & 0 \\ 0 & -3/4 & -1/4 & 1 & -1/4 & 1 \\ 0 & -1/4 & 1/4 & 0 & 1/4 & 0 \\ 0 & -1/4 & 1/4 & 0 & -1/4 & -1/2 \\ 0 & -1/4 & 1/4 & 0 & -1/4 & 1/2 \end{pmatrix}$

Not all of the M matrix is the same for the y and z symmetry Herman-Wallis, but the rows for $D_{10}$ and $D_{11}$ are the same provided that all the cartesian operators have a cyclic permutation of [xyz] their projections.

## SPECIAL CONSIDERATIONS FOR LINEAR MOLECULES

This program set will calculate a variety of interactions and transitions within a Hund's case(b) basis, including spin orbit interactions which change spin multiplicity. The operator $N_a$ in the asymmetric rotor becomes $\lambda$ (or $l$ quantum for bending) for the linear molecule. N is the sum of rotational and electronic orbital angular momenta. For linear molecules, it is convenient (but not essential) to think of the angular momentum along the bond as being purely electronic in nature. In the asymmetric rotor language of this program, the first-order spin orbit interaction takes the operator form of a vector dot product of a direction cosine with the spin vector. It can have two distinct symmetries: $S_a \phi_a = S_a$ connects states of the same lambda, while $S_b \phi_b = S_b$ and $S_c \phi_c = S_c$ connect states where lambda differs by one. When the spin orbit operator connects different spin multiplicity, the reduced matrix value of $< S||S||S' >$ is set to unity.

Use of the symmetries in this program takes some care, particularly for linear molecules where it may not be immediately obvious whether to use the b or the c axis to designate perpendicular operators. For consistency with the parity designation for the symmetric top quanta, the vibronic wave function should be chosen so that it is symmetric with respect to the ab plane. Then the b axis can be used for the inversion-defining axis (i.e. IAX = 2 in the option lines of the .par and .var files can be used to define selection rules under inversion). With this choice, the symmetry of rotation lines in the $D_2$ group is:

|     |                                                              |
| :-- | :----------------------------------------------------------- |
| A   | even N for $\Sigma^-$ states and all N for even $\lambda$ (even parity $\lambda$ doublet) |
| B(a)| odd N for $\Sigma^-$ states and all N for even $\lambda$ (odd parity $\lambda$ doublet)  |
| B(b)| all N for all odd $\lambda$ (even parity $\lambda$ doublet)                  |
| B(c)| all N for all odd $\lambda$ (odd parity $\lambda$ doublet)                   |

For $\Sigma^+$ states, the parity is odd for odd N, while for $\Sigma^-$ states the parity is odd for even N. This means that the Hamiltonian can couple $\Sigma^+$ with $\Sigma^-$ via an operator of B(a) symmetry or B(c) symmetry. An example is an operator like 10200001, which is the $S_a$ spin orbit interaction operator between state v = 0 and v = 1. For normal coupling not involving $\Sigma^-$ states (or for coupling between $\Sigma^-$ states) the operators should have A or B(b) symmetry. Similarly, electric dipole transitions with $\Delta\lambda$ even should have B(a) symmetry, and transitions with $\Delta\lambda = \text{odd}$.

The $g$, $u$ symmetry for an electronic state is for the parity of the wave function under inversion of the space fixed axes. The nuclear exchange symmetry, on the other hand, affects only the statistical weights and does not have any further impact on the factoring of the Hamiltonian. In general, if IAX = 2, WTPL will be the nuclear spin weight for the A and B(b) states, while WTMN will be the weight for the other two symmetries. For $\Sigma_g^+, \Sigma_u^-$, and $g$ states with other $\lambda$, WTPL is the weight for even permutations, while for $\Sigma_u^+, \Sigma_g^-$, and u states with other $\Lambda$, WTPL is the weight for odd permutations. For example, in oxygen, $\Sigma_g^+$ and $\Delta_g$ have WTPL=1 and WTMN=0, while $\Sigma_u^-$ and $\Delta_u$ have WTPL=0 and WTMN=1. The dipole types given above provide for both allowed and forbidden transitions. For transitions that owe their intensity to spin orbit interactions, the effective transition moment will be the product of the interaction and a transition moment to some intermediate state. Examples are a $\Sigma-\Delta$ electric dipole code of 21vv’0 and a magnetic dipole code of 11vv'1. Electric dipole transitions between $\Sigma^+$ and $\Sigma^-$ will use 11vv'0, while magnetic transitions will use 1vv'1. For transitions between $\Sigma^+$ and $\Sigma^+$ the roles of these operators are reversed. Note that for $\Sigma\Sigma$ transitions, 11vv'0 has selection rules of $\Delta N = -2, 0, +2$, while 1vv'1 has selection rules of $\Delta N = -1, +1$.

The correlation between parity and $e$, $f$ designations follow the recommendations of J. M. Brown et al., J. Mol. Spect. **55**, 500 (1975).

|   | odd spin multiplicity | even spin multiplicity |
| :-: | :-------------------- | :--------------------- |
| e | $p = (-1)^{J-1/2}$      | $p = (-1)^{J-1}$       |
| f | $p = (-1)^{J+1/2}$      | $p = (-1)^J$           |


An example of .var file for an oxygen-like molecule:

```
mock oxygen states
4
-3 3 0 0 0 2 0 1 -1 /default for v=0 triplet Sigma-(g)
1 1 2 2 0 2 1 0 -1 /v=1 is singlet Delta(g)
1 2 0 0 0 2 1 0 0 /v=2 is singlet Sigma-(g)
11 1e+7 0 /term value for v=1
22 2e+7 0 /term value for v=2
110010000 1000.0 0 /spin-spin interaction for v=0
99 10000.0 0 /B for all v
```

An example of .int file for an oxygen-like molecule:

```
mock oxygen rotational and electronic transitions
101 99000 200. 0 6 -80. -80. 99999999.
1000 1. /magnetic moment v=0
1110 1. /magnetic moment v=1
11011 1. /sigma - delta magnetic moment
1021 1. /sigma - sigma magnetic moment
```

Note that the forward slash indicates the remainder of the line is a human-readable **comment** and should be ignored by the software.

The quantum number correlations between Hund's case (b) and case (a) can be a bit confusing at first. For example in a doublet $\Pi$ state N=J-1/2 always correlates with $\Omega =3/2$ and N=J+1/2 always correlates with $\Omega =1/2$ on the basis of projection. For A < 0, e.g. OH, the projection-based correlation follows the energy ordering. For A > 0, the lower energy state is $|N=J+1/2, \Omega =3/2\rangle$ as long as $J + 1/2 < \sqrt{|A|/2B}$. Above this J, $|N = J+1/2, \Omega =1/2\rangle$ (based on projection) is higher in energy. Therefore quantum number assignments based on projections lead to different quanta than those based on energy. For a triplet $\Sigma$ state, N=J+1 correlates with $\Sigma = 0$ based on projection, N=J correlates with an odd combination of $\Sigma = 1$ and $\Sigma = -1$, and N=J-1 correlates with an even combination of $\Sigma = 1$ and $\Sigma = -1$.

Since $q$ multiplies the same operator as (B-C) / 2, it is possible to use the sign of $q$ to determine whether there are more electrons in the ab plane ($q > 0$) or whether there are more electrons in the ac plane ($q < 0$).

Explicit approximate relationships for the parameters are:

|           |                  |
| :-------- | :--------------- |
| 100100vv' | A                |
| 100101vv' | $2 A_J$          |
| 1vv'      | B                |
| 100400vv' | $-p/2$           |
| 400vv'    | $q/2$            |
| 200100vv' | a                |
| 1200100vv'| c                |
| 1200000vv'| $(b+c)/3$        |
| 1200400vv'| $-d/2$           |
| 2200100vv'| $1.5 eQq_1$      |
| 2200400vv'| $-eQq_2/4$       |
| 1100100vv'| $4\lambda S(2S-1)$ |

The extra factors of S in the definition of the spin-spin interaction parameter $\lambda$, i.e. a spin-spin interaction $2\lambda (S_z^2 - S^2/3)$, is a correction for the special normalization assumed for $eQq_1$.

## SPECIAL CONSIDERATIONS FOR $l$-doubled STATES

There are three common situations: (1) $l$-doubled states in linear molecules, (2) an n-fold internal rotor with no operators connecting states of differing $\sigma$ and all K having the same nuclear spin states, (3) an n-fold symmetric rotor in a degenerate excited state. In case (1), $l$ correlates with K in the asymmetric rotor and no special basis is needed. See [Special Considerations for Linear Molecules](#special-considerations-for-linear-molecules). For the other two cases, the $l$-doubled states must be specified in adjacent pairs using EWT1. In case (2), the states with values of ($\sigma$ mod n) = 0 or $n/2$ can be treated as $l = 0$ states, and the other degenerate $\sigma$ can be treated as l = $\pm 1$ states. Values of MOD(EWT1, n) = 0 or = n / 2 designate special pairs of states that may be useful for more complicated situations involving two vibronic angular momentum that sum to 0 or n / 2. The states with $l$ = MOD(EWT1, n) < n / 2 (and $l$ > 0) are those with $Kl > 0$. Paired states with EWT1' = n - EWT1 are those with $Kl \le 0$. The sign of K represents the parity, as in the non $l$-doubled states.

For operators that are off-diagonal in $l$, the $\Delta K > 0$ and $\Delta K < 0$ operators are independent. To distinguish these cases, $\Delta K \cdot \Delta l > 0$ if IV2 < IV1 and $\Delta K \cdot \Delta l < 0$ if IV2 > IV1.

For operators with $\Delta l = 0$, two parameter designations for the presence of an $l_z$ operator are possible depending on the state of the NEWLZ flag. If NEWLZ = 1 (preferred), the $l_z$ operator is explicitly specified by designating a parameter with a vibrational state having $l < 0$. The effect of $l_z$ is that the operator for both l-doubled states is multiplied by the sign of l and is specified by a single parameter. Operators with no $l_z$ are specified using the vibrational state with $l > 0$. In this case, the operator is the same for both l-doubled states and is specified by a single parameter. A special example is the 80 series of parameters where 8120vv specifies the average energy of K = 12 l-doubled pair if v is a state with $l = l' > 0$. If v is a state with $l = l' < 0$, the parameter specifies half the difference of the energy of the l-doubled pair. For the 80 series, the behavior is the same when NEWLZ = 0.

If NEWLZ = 0 (the default), the $l_z$ operator is implied by the symmetry of the base operator and is only designated for the state with $l > 0$. Parameters with $l < 0$ are ignored. For a prolate rotor, if EWT2 = EWT2' for two states and the base operator symmetry is a or b, then $l_z$ is present and the symmetry is transformed to null or c, respectively. If EWT2 $\neq$ EWT2' and and the overall symmetry is c or null, then $l_z$ is present and the symmetry is transformed to b or a respectively. For an oblate rotor, a and c are reversed. For both values of NEWLZ, the operator generated is the same. As an example, $P_a P_b$ has c symmetry and would not be multiplied implicitly by $l_z$, but $P_a P_b \sin 2\pi\rho K/3$ has base symmetry of b, is ordinarily imaginary, and would be multiplied by $l_z$. Multiplying by $l_z$, makes the operator have c symmetry and is ordinarily real.

DIAG = 0 is not recommended on the first option line in the .par file, since the first-order energy is not likely to be ordered with $K$.

The $K$ quantum numbers for $l$-doubled states are designated specially when asymmetric rotor quanta are used so that the lower $K$ doublet is associated with the $l > 0$ state and the upper $K$ doublet and $K = 0$ states are associated with $l < 0$. In this way the degenerate states have the same quantum numbers.

With an n-fold symmetric top specified by IAX, operators obey the selection rule that $K-l$ can only change by multiples of n. If there are operators diagonal in $l$ with $\Delta K$ that is not a multiple of n, then the operator is rejected.

## SPECIAL CONSIDERATIONS FOR SYMMETRIC TOPS

The nuclear spin statistics for symmetric top can be specified in a unique way up to 6-fold groups. Start by considering the $D_n$ point groups that have n-fold rotations around the z axis and n $C_2$ rotations perpendicular to the z axis. For bosons and for even sets of equivalent fermions, the symmetry of the overall wave function must have $A_1$ symmetry. For one set of equivalent fermions, the overall symmetry is $A_2, B_2, A_1, B_1$ for n = 3, 4, 5, 6, respectively. The influence of vibronic symmetry is based on user input. $B_1$ and $B_2$ symmetries are specified in separate option lines with ABS(IAX) = 5, 8, 11. The statistical weights for each vibronic state are specified by WTPL, WTMN, and EWT0. For even n, the user specifies two option lines with different values of IAX in order to cover all symmetries. For $D_4$ symmetry EWT0 is ignored for IAX = 8. For ABS(IAX) > 3, unspecified weights are assumed to be zero. The symmetry of the vibronic function is defined by $l$ and the symmetry under the perpendicular $C_2$ operator. The symmetry of the rotational wave functions is $A_1$ or $A_2$ for MOD($K-l$, n) = 0, $B_1$ or $B_2$ for MOD($K-l$, n) = n/2, and E otherwise. (Two-fold rotations are handled in a straight-forward way using the full symmetry of the $D_2$ group and will not be considered further in this section.) The two-fold symmetry (e.g. $A_1$ vs. $A_2$) is defined by the Wang symmetry along the b axis and atom 1 is assumed to pass through this perpendicular $C_2$ axis. The nuclear spin statistics are readily defined in terms of the permutation symmetry of the nuclei. The values given below are for n equivalent nuclei with spin I that lie in a plane perpendicular to the z axis. N is the number of nuclear spin states with a given symmetry. The number following $\Rightarrow$ is the value for $I = 1/2$.

For $D_2$ symmetry:

- $N(A_1) = (2I+1)(I+1) \Rightarrow 3$
- $N(B_x) = 0$
- $N(B_y) = 0$
- $N(B_z) = (2I+1)(I) \Rightarrow 1$

For $D_3$ symmetry:

- $N(A_1) = (2I+1) [(2I+1)^2 + 3(2I+1)+2]/6 \Rightarrow 4$
- $N(A_2) = (2I+1) [(2I+1)^2 - 3(2I+1)+2]/6 \Rightarrow 0$
- $N(E) = (2I+1)(I+1)(8I)/3 \Rightarrow 4$

For $D_4$ symmetry:

- $N(A_1) = (2I+1) [(2I+1)^3 + 2(2I+1)^2 + 3(2I+1) + 2]/8 \Rightarrow 6$
- $N(A_2) = (2I+1) [(2I+1)^3 - 2(2I+1)^2 - (2I+1) + 2]/8 \Rightarrow 0$
- $N(B_1) = (2I+1) [(2I+1)^3 - 2(2I+1)^2 + 3(2I+1) - 2]/8 \Rightarrow 1$
- $N(B_2) = (2I+1) [(2I+1)^3 + 2(2I+1)^2 - (2I+1) - 2]/8 \Rightarrow 3$
- $N(E) = (2I+1)^2(I+1)(2I) \Rightarrow 6$

For $D_5$ symmetry:

- $N(A_1) = (2I+1) [(2I+1)^4 + 5(2I+1)^2 + 4]/10 \Rightarrow 8$
- $N(A_2) = (2I+1) [(2I+1)^4 - 5(2I+1)^2 + 4]/10 \Rightarrow 0$
- $N(E_1) = (2I+1) [(2I+1)^4 - 1]2/5 \Rightarrow 12$
- $N(E_2) = N(E_1)$

For $D_6$ symmetry:

- $N(A_1) = (2I+1) [(2I+1)^5 + 3(2I+1)^3 + 4(2I+1)^2 + 2(2I+1) + 2]/12 \Rightarrow 13$
- $N(A_2) = (2I+1) [(2I+1)^5 - 3(2I+1)^3 - 2(2I+1)^2 + 2(2I+1) + 2]/12 \Rightarrow 1$
- $N(B_1) = (2I+1) [(2I+1)^5 + 3(2I+1)^3 - 4(2I+1)^2 + 2(2I+1) - 2]/12 \Rightarrow 7$
- $N(B_2) = (2I+1) [(2I+1)^5 - 3(2I+1)^3 + 2(2I+1)^2 + 2(2I+1) - 2]/12 \Rightarrow 3$
- $N(E_1) = (2I+1) [(2I+1)^5 - (2I+1)^2 - (2I+1) + 1]/3 \Rightarrow 18$
- $N(E_2) = (2I+1) [(2I+1)^5 + (2I+1)^2 - (2I+1) - 1]/3 \Rightarrow 22$

For $C_n$ symmetry: $N(A)' = N(A_1) + N(A_2)$, $N(B)' = N(B_1) + N(B_2)$, and $N(E)' = 2 N(E)$. Nuclear spin statistics of symmetries involving the inversion operator, reflection operator, or an improper rotation have the same nuclear spin statistics as the pure rotation sub-group. The symmetry of the rotational state being considered is then multiplied by the nuclear spin symmetry that gives the correct overall symmetry.

For $D_n$ symmetry, WTPL = $N(A_1)$ or $N(B_1)$ for bosons in $A_1$ vibronic states and WTPL = $N(A_2)$ or $N(B_2)$ for fermions in $A_1$ vibronic states. WTMN = $N(A_2)$ or $N(B_2)$ for bosons in $A_1$ vibronic states and WTMN = $N(A_1)$ or $N(B_1)$ for fermions in $A_1$ vibronic states. EWT0 is the weight for E rotational states ($E_1$ or $E_2$ for $D_6$ symmetry). Due to the degeneracy of the spin states, EWT0 = N(E) / 4. For the special case of n = 3 when WTPL = 0 or WTMN = 0, EWT0 should be doubled because only half of the states with MOD(K-l) $\neq$ 0 are calculated. For non-degenerate vibronic states, values of WTPL and WTMN may need to be reversed. For a E vibrational state, use of l-doubling allows the same statistical weights as the corresponding non-degenerate vibrational states. For $C_n$ symmetry, WTPL = WTMN = N(A)' or N(B)' and EWT0 = N(E) / 4.

As an example consider $NH_3$ with $D_3$ rotational symmetry and overall symmetry of $A_2$. For an $A_1$ vibrational state, WTPL = $N(A_2)$ = 0, WTMN = $N(A_1)$ = 4, and EWT0 = 2 N(E) / 4 = 2. For an $A_2$ vibrational state, WTPL and WTMN are reversed. For $ND_3$ with $D_3$ rotational symmetry, the overall symmetry is $A_1$. For an $A_1$ vibrational state, WTPL = $N(A_1)$ = 10, WTMN = $N(A_2)$ = 1, and EWT0 = N(E) / 4 = 4. For $CH_3F$ with $C_3$ rotational symmetry and overall symmetry of A, WTPL = WTMN = N(A)' = 4 and EWT0 = N(E)' / 4 = 2.

For reference, the symmetry multiplication table is:

|           |           |           |           |               |               |
| :-------- | :-------- | :-------- | :-------- | :-----------: | :-----------: |
| $A_1$     | $A_2$     | $B_1$     | $B_2$     | $E_1$         | $E_2$         |
| $A_2$     | $A_1$     | $B_2$     | $B_1$     | $E_1$         | $E_2$         |
| $B_1$     | $B_2$     | $A_1$     | $A_2$     | $E_2$         | $E_1$         |
| $B_2$     | $B_1$     | $A_2$     | $A_1$     | $E_2$         | $E_1$         |
| $E_1$     | $E_1$     | $E_2$     | $E_2$     | $A_1+A_2+E_2$ | $B_1+B_2+E_1$ |
| $E_2$     | $E_2$     | $E_1$     | $E_1$     | $B_1+B_2+E_1$ | $A_1+A_2+E_2$ |

For n = 3, E $\times$ E = $A_1 + A_2 + E$. For n = 4, E $\times$ E = $A_1 + A_2 + B_1 + B_2$. Whenever only one resultant of a E $\times$ E product is allowed, N(E) must be divided by 4.

For hyperfine interactions of the n nuclei in an n-fold top, the user can either use standard coupling or (preferably) the $I_{tot}$ basis indicated by the sign of IAX. The standard coupling is more flexible, but the dark states with inappropriate overall symmetry are not readily identifiable. In the $I_{tot}$ basis, only parameters with $I_1$ or $[I_1 \times I_2]$ need be specified. The actual spin operator used is $I_a = I_1 + \epsilon I_2 + \epsilon^2 I_3 + \ldots$, where $\epsilon = \text{EXP}(2\pi i/n)$. The value of $\alpha$ is given implicitly by MOD($\Delta I - \Delta K, n$). A similar Fourier series is provided for $[I_1 \times I_2]$ (for L = 0,2) and for quadrupole interactions, namely $[I_1 \times I_1]^{(2)}$. Note that $\alpha = 0$ operators couple spin states with the same symmetry, while other values of $\alpha$ couple ortho and para. For $D_n$ symmetry and the $I_{tot}$ basis, EWT0 = 0 and either WTPL or WTMN is zero based on the overall symmetry required. For $C_n$ symmetry and the $I_{tot}$ basis, EWT0 = 0 and WTPL and WTMN are both zero or one based on the overall symmetry required. The matrix elements of the operators are constructed so that there is no connection between states of different overall symmetry.

## Asymmetric rotor example

Example of parameter types for asymmetric rotors (assuming < 10 vibronic states):

|          |                                                   |
| -------: | :------------------------------------------------ |
| 11       | energy for v=1                                    |
| 10000    | $A_{00}$                                          |
| 10099    | A (all vibrational states with v'=v")             |
| 20099    | B                                                 |
| 30099    | C                                                 |
| 40099    | $0.25*(B-C)$ (if prolate basis selected)          |
| 299      | $-D_J$                                            |
| 1199     | $-D_{JK}$                                         |
| 2000     | $-D_{K00}$                                        |
| 600001   | i $N_z$ interaction between v=0 and v=1           |
| 20000099 | NI for second spin                                |
| 120010099| $S_a I_a$                                         |
| 220010099| $1.5 *$ quadrupole moment $c_{zz}$ for second spin |
| 220040099| Quadrupole moment $0.25*(c_{xx}-c_{yy})$ for second spin |

Quadrupole and magnetic spin-spin interactions are defined to be traceless (i.e. $c_{xx} + c_{yy} + c_{zz} = 0$ or $T_{xx} + T_{yy} + T_{zz} = 0$) and all three components cannot be simultaneously fit. The most efficient choice of parameters is shown in the table above. In cases where the user wants an alternative, it is possible to use constrained parameters. For example, to fit $c_{aa}$ and $c_{cc}$ (with no multipliers):

```
220010099 100.
-220020099 -100.
220030099 50.
-220020099 -50.
```

specifies $\chi_{aa} = 100$, $\chi_{cc} = 50$, and $\chi_{bb} = -150$.
