#ifndef NONDIMENSIONALNUMBERS_H
#define NONDIMENSIONALNUMBERS_H

/* Definition of non-dimensional numbers */
double Re; // Reynolds Number
double Sc; // schmidt Number
double Nu; // Nusselt Number
double Le; // Lewis Number
double BM; // Spalding mass transfer number
enum spaldingComputationMethod
{
    massFraction,
    moleFraction
};
spaldingComputationMethod spaldingMethod;

double Gr; // Grashof Number
double Pr; // Prandtl Number
double Ra; // Rayleigh Number
double Sh; // Sherwood Number
enum convectionComputationMethod
{
    naturalConvection,
    turbulentForcedConvection
};
convectionComputationMethod convectionMethod;

#endif