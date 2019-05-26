# STCSIR
Repository for the stable time-changed SIR model

Description of the files:
MLgen.R contains the function for the generation of a vector of Mittag-Leffler random variables;
SIRsimexp.R contains the function for the simulation of a classical SIR model;
SIRsimfrac.R contains the function for the simulation of a time-changed SIR model;
Figure1.R and Figure2.R are some applications of these functions.

To use the functions, the package "stabledist" is needed. Moreover, to use the functions in SIRsimexp.R and SIRsimfrac.R, it is necessary to import first MLgen.R.
