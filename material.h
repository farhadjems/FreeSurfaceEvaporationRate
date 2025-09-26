#ifndef MATERIAL_H
#define MATERIAL_H

#include "constant.h"
/* Definition of variable properties coefficients*/
struct VariableProperties
{
    double k;
    double Cp;
    double mu;
};

// Fluid

enum FluidName // The type of fluid in simulation
{
    water,
    heptane,
    acetone,
    methanol,
    userFluid,
};

struct Fluid
{
    FluidName name;
    double Tb; // Boiling Temperature of fluid - unit: Kelvin
    double Lv; // Latent heat of boiling - unit: kJ/kg
    double MW; // Molar weight - unit: kg/kmol
    // Liquid phase Properties
    double kLiq;    // Thermal conductivity - unit: kW/(m.K)
    double CpLiq;   // Thermal specific heat capacity of constant pressure - unit: kJ/(kg.K)
    double muLiq;   // Viscosity - unit: Pa.s
    double rhoLiq;  // Density - unit: kg/m3
    double nuLiq;   // Kinematic Viscosity - unit: m3/s2
    double eL;      // Emissivity - non-dimensional
    double alpha;   // Absorptivity - non-dimensional
    double mLiq;    // Mass of liquid (fuel) - unit: kg
    double evpRate; // Evaporation rate of liquid in free surface per unit area (MLRPUA) - unit: kg/s
    double dLiq;    // Liquid depth - unit: m
    // Vapor Properties
    double kVap;             // Thermal conductivity - unit: kW/(m.K)
    double CpVap;            // Thermal specific heat capacity of constant pressure - unit: kJ/(kg.K)
    double muVap;            // Viscosity - unit: Pa.s
    double rhoVap;           // Density - unit: kg/m3
    double nuVap;            // Kinematic Viscosity - unit: m3/s2
    double Ysurf;            // Mass fraction of vapor at surface of liquid pool
    double Xsurf;            // Mole fraction of vapor at surface of liquid pool
    double Yfilm;            // Mass fraction of vapor at film layer
    double Xfilm;            // Mole fraction of vapor at film layer
    double Ysurr;            // Mass fraction of vapor at surrounding
    double Xsurr;            // Mole fraction of vapor at surrounding
    VariableProperties A[6]; // coefficients of variable properties
};

// Gas
struct Gas
{
    double MW;               // Molar weight - unit: kg/kmol
    double k;                // Thermal conductivity - unit: kW/(m.K)
    double Cp;               // Thermal specific heat capacity of constant pressure - unit: kj/(kg.K)
    double mu;               // Viscosity - unit: Pa.s
    double rho;              // Density - unit: kg/m3
    double nu;               // Kinematic Viscosity - unit: m3/s2
    double Ysurf;            // Mass fraction of vapor at surface of liquid pool
    double Xsurf;            // Mole fraction of vapor at surface of liquid pool
    double Yfilm;            // Mass fraction of vapor at film layer
    double Xfilm;            // Mole fraction of vapor at film layer
    double Ysurr;            // Mass fraction of vapor at surrounding
    double Xsurr;            // Mole fraction of vapor at surrounding
    VariableProperties A[6]; // coefficients of variable properties
};

// Defining Functions
void gasVariableProperties(Gas &gas, double TK);
void gasVariablePropertiesCoefficients(Gas &gas, double TK);
void gasConstantProperties(Gas &gas);
void fluidVariableProperties(Fluid &fluid, double TK);
void fluidVariablePropertiesCoefficients(Fluid &fluid, double TK);
void fluidConstantProperties(Fluid &fluid);

#endif