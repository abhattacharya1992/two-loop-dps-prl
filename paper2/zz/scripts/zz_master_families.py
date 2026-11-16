"""Identify the finite FORM export's labels in the primary master basis.

The suffix c exchanges t and u. The suffix m exchanges ma2 and mb2 while
keeping t fixed, so m is x12*x34 and cm is x34. In particular, MN13m and
MN13cm must not be interpreted as N13x34 and N13x34x12 respectively.
The nonplanar weight-four coefficients of the unrenormalized export provide
an independent check of these identifications, without using C++ amplitudes.
"""

FAMILY = {
    "MP12": "P12", "MP12c": "P12x12",
    "MP13": "P13", "MP13c": "P13x12",
    "MP23": "P23", "MP23c": "P23x12",
    "Maddon": "addon", "Maddonc": "addonx12",
    "MN12": "N12", "MN13": "N13", "MN13c": "N13x12",
    "MN13m": "N13x34x12", "MN13cm": "N13x34", "MN34": "N34",
}
