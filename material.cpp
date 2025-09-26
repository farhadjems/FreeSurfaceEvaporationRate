#include "material.h"

/* Definition of Property functions*/

// Fluid Properties
void fluidConstantProperties(Fluid &fluid)
{
    // Table 1 - Liquid Properties
    switch (fluid.name)
    {
    case water:                // H2O
        fluid.MW = 18.0;       // unit: kg/kmol
        fluid.rhoLiq = 1000.0; // unit: kg/m3
        fluid.CpLiq = 4.18;    // unit: kj/(kg.k)
        fluid.Lv = 2260.0;     // unit: kj/kg
        fluid.eL = 0.95;       // non-dimensional
        fluid.Tb = 373.0;      // unit: Kelvin
        fluid.alpha = 1.0;     // non-dimensional
        break;

    case heptane:             // C7H16
        fluid.MW = 100.0;     // unit: kg/kmol
        fluid.rhoLiq = 675.0; // unit: kg/m3
        fluid.CpLiq = 2.24;   // unit: kj/(kg.k)
        fluid.Lv = 317.0;     // unit: kj/kg
        fluid.eL = 0.95;      // non-dimensional
        fluid.Tb = 371.5;     // unit: Kelvin
        fluid.alpha = 1.0;    // non-dimensional
        break;

    case acetone:             // C3H6O
        fluid.MW = 58.0;      // unit: kg/kmol
        fluid.rhoLiq = 791.0; // unit: kg/m3
        fluid.CpLiq = 2.13;   // unit: kj/(kg.k)
        fluid.Lv = 501.0;     // unit: kj/kg
        fluid.eL = 0.95;      // non-dimensional
        fluid.Tb = 329.3;     // unit: Kelvin
        fluid.alpha = 1.0;    // non-dimensional
        break;

    case methanol:            // CH3OH
        fluid.MW = 32.0;      // unit: kg/kmol
        fluid.rhoLiq = 796.0; // unit: kg/m3
        fluid.CpLiq = 2.48;   // unit: kj/(kg.k)
        fluid.Lv = 1099.0;    // unit: kj/kg
        fluid.eL = 0.95;      // non-dimensional
        fluid.Tb = 337.8;     // unit: Kelvin
        fluid.alpha = 1.0;    // non-dimensional
        break;

    case userFluid:         // user-defined material. enter the properties and change the variable "liquidName" value to "user"
        fluid.MW = 1.0;     // unit: kg/kmol
        fluid.rhoLiq = 1.0; // unit: kg/m3
        fluid.CpLiq = 1.0;  // unit: kj/(kg.k)
        fluid.Lv = 1.0;     // unit: kj/kg
        fluid.eL = 1.0;     // non-dimensional
        fluid.Tb = 273.0;   // unit: Kelvin
        fluid.alpha = 1.0;  // non-dimensional
        break;
    }
}

void fluidVariablePropertiesCoefficients(Fluid &fluid, double TK)
{
    double TC = TK - 273.0;
    switch (fluid.name)
    {
    case water:
        // vap.Cp coefficients
        fluid.A[0].Cp = 2.67703787;
        fluid.A[1].Cp = 2.97318329e-3;
        fluid.A[2].Cp = -7.73769690e-7;
        fluid.A[3].Cp = 9.44336689e-11;
        fluid.A[4].Cp = -4.26900959e-15;
        fluid.A[5].Cp = 0;
        // vap.k and vap.mu coefficients
        if (TK >= 295 && TK <= 1000)
        {
            fluid.A[0].k = 0.15541443e1;
            fluid.A[1].k = 0.66106305e2;
            fluid.A[2].k = 0.55968860e4;
            fluid.A[3].k = -0.39259598e1;
            fluid.A[4].k = 0.0;
            fluid.A[5].k = 0.0;

            fluid.A[0].mu = 0.78387780e0;
            fluid.A[1].mu = -0.38260408e3;
            fluid.A[2].mu = 0.49040158e5;
            fluid.A[3].mu = 0.85222785e0;
            fluid.A[4].mu = 0.0;
            fluid.A[5].mu = 0.0;
        }
        if (TK > 1000 && TK <= 5000)
        {
            fluid.A[0].k = 0.79349503e0;
            fluid.A[1].k = -0.13340063e4;
            fluid.A[2].k = 0.378664327e6;
            fluid.A[3].k = 0.23591474e1;
            fluid.A[4].k = 0.0;
            fluid.A[5].k = 0.0;

            fluid.A[0].mu = 0.50714993e0;
            fluid.A[1].mu = -0.68966913e3;
            fluid.A[2].mu = 0.874547505e5;
            fluid.A[3].mu = 0.30285155e1;
            fluid.A[4].mu = 0.0;
            fluid.A[5].mu = 0.0;
        }
        break;

    case heptane:
        // vap.Cp coefficients
        fluid.A[0].Cp = -1.3934e-1;
        fluid.A[1].Cp = 2.6295e-3;
        fluid.A[2].Cp = 1.15e-5;
        fluid.A[3].Cp = -1.7062e-7;
        fluid.A[4].Cp = -1.205e-9;
        fluid.A[5].Cp = 1.5923e-11;
        // vap.k coefficients
        fluid.A[0].k = -4.6160e0;
        fluid.A[1].k = 7.7908e-3;
        fluid.A[2].k = -1.5605e-5;
        fluid.A[3].k = 1.158e-7;
        fluid.A[4].k = -1.4660e-9;
        fluid.A[5].k = 5.4467e-12;
        // vap.mu coefficients
        fluid.A[0].mu = 3.8877;
        fluid.A[1].mu = 6.1806e-3;
        fluid.A[2].mu = -3.0481e-5;
        fluid.A[3].mu = 9.7703e-8;
        fluid.A[4].mu = 0.0;
        fluid.A[5].mu = 0.0;
        break;

    case acetone:
        // vap.Cp coefficients
        fluid.A[0].Cp = 1.9452e-1;
        fluid.A[1].Cp = 2.2900e-3;
        fluid.A[2].Cp = -8.6330e-7;
        fluid.A[3].Cp = -2.0672e-8;
        fluid.A[4].Cp = 1.9250e-10;
        fluid.A[5].Cp = -6.5969e-13;
        // vap.k coefficients
        fluid.A[0].k = -4.6467;
        fluid.A[1].k = 7.1871e-3;
        fluid.A[2].k = -2.0976e-5;
        fluid.A[3].k = 4.6070e-7;
        fluid.A[4].k = -5.4286e-9;
        fluid.A[5].k = 1.9213e-11;
        // vap.mu coefficients
        fluid.A[0].mu = 4.3519e0;
        fluid.A[1].mu = 2.6016e-3;
        fluid.A[2].mu = -3.5279e-6;
        fluid.A[3].mu = 2.1574e-7;
        fluid.A[4].mu = -3.5172e-9;
        fluid.A[5].mu = 1.3856e-11;
        break;

    case methanol:
        // vap.Cp coefficients
        fluid.A[0].Cp = 0.31404e0;
        fluid.A[1].Cp = 5.01e-3;
        fluid.A[2].Cp = -1.7365e-4;
        fluid.A[3].Cp = 3.4375e-6;
        fluid.A[4].Cp = -2.5001e-8;
        fluid.A[5].Cp = 5.9705e-11;
        // vap.k coefficients
        fluid.A[0].k = -6.7646e0;
        fluid.A[1].k = 8.9038e-3;
        fluid.A[2].k = -3.5238e-5;
        fluid.A[3].k = 1.4232e-7;
        fluid.A[4].k = -4.003e-10;
        fluid.A[5].k = 5.0729e-13;
        // vap.mu coefficients
        fluid.A[0].mu = 4.4794e0;
        fluid.A[1].mu = 3.3930e-3;
        fluid.A[2].mu = 4.0163e-6;
        fluid.A[3].mu = -7.3376e-8;
        fluid.A[4].mu = 2.4641e-10;
        fluid.A[5].mu = -2.1042e-13;
        break;

    case userFluid: // Enter the constants here
        // vap.Cp coefficients
        fluid.A[0].Cp = 1.0;
        fluid.A[1].Cp = 1.0;
        fluid.A[2].Cp = 1.0;
        fluid.A[3].Cp = 1.0;
        fluid.A[4].Cp = 1.0;
        fluid.A[5].Cp = 1.0;
        // vap.k coefficients
        fluid.A[0].k = 1.0;
        fluid.A[1].k = 1.0;
        fluid.A[2].k = 1.0;
        fluid.A[3].k = 1.0;
        fluid.A[4].k = 1.0;
        fluid.A[5].k = 1.0;
        // vap.mu coefficients
        fluid.A[0].mu = 1.0;
        fluid.A[1].mu = 1.0;
        fluid.A[2].mu = 1.0;
        fluid.A[3].mu = 1.0;
        fluid.A[4].mu = 1.0;
        fluid.A[5].mu = 1.0;
        break;
    }
}

void fluidVariableProperties(Fluid &fluid, double TK)
{
    double TC = TK - 273.0;
    fluidVariablePropertiesCoefficients(fluid, TK);
    switch (fluid.name)
    {
    case water:
        fluid.CpVap = (fluid.A[0].Cp + fluid.A[1].Cp * TK + fluid.A[2].Cp * pow(TK, 2) + fluid.A[3].Cp * pow(TK, 3) + fluid.A[4].Cp * pow(TK, 4)) * R / fluid.MW; // Equation A1
        fluid.kVap = exp(fluid.A[0].k * log(TK) + fluid.A[1].k / TK + fluid.A[2].k / pow(TK, 2) + fluid.A[3].k) * 1e-7;                                           // Equation A2
        fluid.muVap = exp(fluid.A[0].mu * log(TK) + fluid.A[1].mu / TK + fluid.A[2].mu / pow(TK, 2) + fluid.A[3].mu) * 1e-7;                                      // Equation A3
        break;

    case (heptane):
        fluid.CpVap = exp(fluid.A[0].Cp + fluid.A[1].Cp * TC + fluid.A[2].Cp * pow(TC, 2) + fluid.A[3].Cp * pow(TC, 3) +
                          fluid.A[4].Cp * pow(TC, 4) + fluid.A[5].Cp * pow(TC, 5));                                                                                                       // Equation A3
        fluid.kVap = exp(fluid.A[0].k + fluid.A[1].k * TC + fluid.A[2].k * pow(TC, 2) + fluid.A[3].k * pow(TC, 3) + fluid.A[4].k * pow(TC, 4) + fluid.A[5].k * pow(TC, 5)) * 1e-3;        // Equation A3
        fluid.muVap = exp(fluid.A[0].mu + fluid.A[1].mu * TC + fluid.A[2].mu * pow(TC, 2) + fluid.A[3].mu * pow(TC, 3) + fluid.A[4].mu * pow(TC, 4) + fluid.A[5].mu * pow(TC, 5)) * 1e-7; // Equation A3
        break;

    case (acetone):
        fluid.CpVap = exp(fluid.A[0].Cp + fluid.A[1].Cp * TC + fluid.A[2].Cp * pow(TC, 2) + fluid.A[3].Cp * pow(TC, 3) + fluid.A[4].Cp * pow(TC, 4) + fluid.A[5].Cp * pow(TC, 5));        // Equation A3
        fluid.kVap = exp(fluid.A[0].k + fluid.A[1].k * TC + fluid.A[2].k * pow(TC, 2) + fluid.A[3].k * pow(TC, 3) + fluid.A[4].k * pow(TC, 4) + fluid.A[5].k * pow(TC, 5)) * 1e-3;        // Equation A3
        fluid.muVap = exp(fluid.A[0].mu + fluid.A[1].mu * TC + fluid.A[2].mu * pow(TC, 2) + fluid.A[3].mu * pow(TC, 3) + fluid.A[4].mu * pow(TC, 4) + fluid.A[5].mu * pow(TC, 5)) * 1e-7; // Equation A3
        break;

    case (methanol):
        fluid.CpVap = exp(fluid.A[0].Cp + fluid.A[1].Cp * TC + fluid.A[2].Cp * pow(TC, 2) + fluid.A[3].Cp * pow(TC, 3) + fluid.A[4].Cp * pow(TC, 4) + fluid.A[5].Cp * pow(TC, 5));        // Equation A3
        fluid.kVap = exp(fluid.A[0].k + fluid.A[1].k * TC + fluid.A[2].k * pow(TC, 2) + fluid.A[3].k * pow(TC, 3) + fluid.A[4].k * pow(TC, 4) + fluid.A[5].k * pow(TC, 5)) * 1e-3;        // Equation A3
        fluid.muVap = exp(fluid.A[0].mu + fluid.A[1].mu * TC + fluid.A[2].mu * pow(TC, 2) + fluid.A[3].mu * pow(TC, 3) + fluid.A[4].mu * pow(TC, 4) + fluid.A[5].mu * pow(TC, 5)) * 1e-7; // Equation A3
        break;

    case userFluid:
        fluid.CpVap = 1.0;
        fluid.kVap = 1.0;
        fluid.muVap = 1.0;
        break;
    }
}

// Air Properties
void gasConstantProperties(Gas &gas)
{
    gas.MW = 29.0; // unit: kg/kmol
}

void gasVariablePropertiesCoefficients(Gas &gas, double TK)
{
    double TC = TK - 273.0;
    // Cp
    gas.A[0].Cp = 7.5512e-3;
    gas.A[1].Cp = 9.8834e-6;
    gas.A[2].Cp = 5.0936e-7;
    gas.A[3].Cp = -5.3667e-10;
    gas.A[4].Cp = 2.1073e-13;
    gas.A[5].Cp = -2.3807e-17;
    // k
    gas.A[0].k = 3.1349;
    gas.A[1].k = 3.5533e-3;
    gas.A[2].k = -5.2910e-6;
    gas.A[3].k = 4.4429e-9;
    gas.A[4].k = -1.7077e-12;
    gas.A[5].k = 2.4576e-16;
    // mu
    gas.A[0].mu = 5.0907;
    gas.A[1].mu = 3.2040e-3;
    gas.A[2].mu = -4.7190e-6;
    gas.A[3].mu = 3.8481e-9;
    gas.A[4].mu = -1.4484e-12;
    gas.A[5].mu = 2.0023e-16;
}

void gasVariableProperties(Gas &gas, double TK)
{
    double TC = TK - 273.0; // Input Temperature - Celsius
    gasVariablePropertiesCoefficients(gas, TK);
    gas.Cp = exp(gas.A[0].Cp + gas.A[1].Cp * TC + gas.A[2].Cp * pow(TC, 2) + gas.A[3].Cp * pow(TC, 3) + gas.A[4].Cp * pow(TC, 4) + gas.A[5].Cp * pow(TC, 5));        // Equation A3
    gas.k = exp(gas.A[0].k + gas.A[1].k * TC + gas.A[2].k * pow(TC, 2) + gas.A[3].k * pow(TC, 3) + gas.A[4].k * pow(TC, 4) + gas.A[5].k * pow(TC, 5)) * 1e-6;        // Equation A3
    gas.mu = exp(gas.A[0].mu + gas.A[1].mu * TC + gas.A[2].mu * pow(TC, 2) + gas.A[3].mu * pow(TC, 3) + gas.A[4].mu * pow(TC, 4) + gas.A[5].mu * pow(TC, 5)) * 1e-7; // Equation A3
    gas.rho = (Patm * gas.MW) / (R * TK);
}
