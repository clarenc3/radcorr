# Radiative corrections
Simple program to calculate size of radiative corrections calculation using input from Oleksandr Tomalak et al. Phys.Rev.D 106 (2022) 9, 093006, also on 
[hep-ph arXiv:2204.11379](https://arxiv.org/abs/2204.11379)

# Dependencies
Uses calculations from `ROOT` file which contains `TGraph`, hence dependent on the `ROOT` framework.

# Compiling
Set up `root`, type `make`. e.g. 
```
source /path/to/root/bin/thisroot.sh
make
```
which should produce a shared library `libRadCorr.so`, an object `RadCorrCalculator.o`, an executable `Calculate.exe`, and an executable `CheckQ2.exe`.

# Testing and validating
The `Calculate.exe` program produces a `pdf` and `root` file showing weights for muon neutrinos and anti-neutrinos.

The `CheckQ2.exe` program produces text to screen with the weights at the Q2 boundary. This can be used for validation:
```
Enu:   0.2 | Q2max:  0.09226 | Weight:   1.0282
Enu:   0.3 | Q2max:  0.20129 | Weight:  0.42489
Enu:   0.5 | Q2max:  0.46755 | Weight:  0.45413
Enu:   0.6 | Q2max:  0.61596 | Weight:        0
Enu:  0.75 | Q2max:  0.85067 | Weight:  0.54474
Enu:   0.9 | Q2max:   1.0957 | Weight:        0
Enu:     1 | Q2max:   1.2631 | Weight:        0
Enu:     2 | Q2max:   3.0284 | Weight:   0.5604
Enu:     3 | Q2max:   4.8586 | Weight:   0.7618
Enu:     5 | Q2max:   8.5712 | Weight:   2.0608
Enu:    10 | Q2max:   17.925 | Weight:   9.8641
Enu:    15 | Q2max:   27.301 | Weight:   5.2827
Enu:    20 | Q2max:   36.684 | Weight:   7.3449
Enu:    25 | Q2max:   46.069 | Weight:   9.5228
Enu:    30 | Q2max:   55.456 | Weight:       10
Enu:    35 | Q2max:   64.843 | Weight:       10
Enu:    40 | Q2max:   74.231 | Weight:       10
Enu:    45 | Q2max:   83.619 | Weight:       10
Enu:    50 | Q2max:   93.007 | Weight:       10
Enu:    70 | Q2max:   130.56 | Weight:       10
```

# Caveats
* Linear interpolation is intended method. Spline interpolation is supported but not validated and used. Expect interesting results, especially around the Q2 boundary.
* The maximum supported neutrino energy is 70 GeV. Anything above that will apply the same weight as 70 GeV.
* The minimum supported neutrino energy is 0.2 GeV. Anything below that will apply the same weight as 0.2 GeV.
* Calculations are not expected to be reliable at higher Q2, e.g. Q2 > 3 GeV2. The inputs are not consistently provided about 3 GeV2, for example the 70 GeV pre-calculated graph stops at 3 GeV2. Above this point in Q2, a linear or spline interpolation will be used, which can give radically different results. If you need calculations above this Q2, please contact Oleksandr Tomalak. An example `root` application is provided to check the impact of the extrpolation, see `CheckExtrapolation.cpp`

# Examples
![numu](https://github.com/clarenc3/radcorr/assets/11278399/4c765c9a-5740-4f18-85ce-bbe80919d8f5) 

![numubar](https://github.com/clarenc3/radcorr/assets/11278399/75170063-b042-45fb-87b2-94c3feea87e5)
