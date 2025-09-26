#ifndef DOMAIN_H
#define DOMAIN_H

/* Definition of domains*/
struct Domain
{
    double MW;  // Molar weight - unit: kg/kmol
    double T;   // Temperature - unit: Kelvin
    double k;   // Thermal conductivity - unit: kW/(m.K)
    double Cp;  // Thermal specific heat capacity of constant pressure - unit: kJ/(kg.K)
    double mu;  // Viscosity - unit: Pa.s
    double rho; // Density - unit: kg/m3
    double nu;  // Kinematic Viscosity - unit: m3/s2
};

/* Definition of Liquid Pool */
struct Pool
{
    double A; // Area - unit: m2
    double P; // Perimeter - unit: m
    double D; // Diameter - unit: m
    double V; // Volume - unit: m3
    double L; // Ratio of area to perimeter - unit: m
};

#endif