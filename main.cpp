#include <ctime>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <sys/stat.h>

#include "material.h"
#include "domain.h"
#include "nonDimensionalNumbers.h"

#define VERSION 3.0
// #define DEBUG

using namespace std;

/********** DEFINITION OF VARIABLES **********/
Fluid fluid;
Gas gas;

Domain surf;
Domain surr;
Domain film;

Pool pool;

// Heat flux per unit area

string heatSource = "";
string SurroundingT = "";
string address = "";
double q; // Heat flux - unit: kW/m2

// Temporal variables
double delt; // Time-step size - unit: seconds
int kt;      // Number of time steps
double t;    // Total simulation time - unit: seconds

// Constants - Equations 35 - 43
double a1;       // Equation 39 - unit: K-3
double a2;       // Equation 40 - unit: K
double a3;       // Equation 35 - unit: non-dimensional
double a4;       // Equation 36 - unit: non-dimensional
double Tamb;     // Ambient temperature in Radiation calculations - unit: Kelvin
double Tstar;    // Equation 37 - unit: Kelvin
double Taw;      // Equation 38 - unit: second
double ddt;      // Equations 42 & 43 - unit: (Kelvin)-1
double surfTnew; // Temperature at new time-step - unit: Kelvin

/********** DEFINITION OF FUNCTIONS **********/
filesystem::path mainPath = filesystem::current_path();
string projectAddress = static_cast<string>(mainPath);
string variableNames = "";

bool convertToEnum(string fluidName);
bool checkInput();
void createFile(string &address);
void saveData(int k, string &address);
void generalSetup();
void initialConditions();
void computeInlets(int k);
void liquidSurface();
void filmProperties();
void nonDimensionalNumbers();
void constants();
void temperature();
void evaporationRate();
void fuelNewValues();
void printData(int k);

/* Main function of code. DO NOT CHANGE IT!!! */
int main(int argc, char *argv[])
{
    heatSource = static_cast<string>(argv[1]);
    SurroundingT = static_cast<string>(argv[2]);
    if (!checkInput())
    {
        return 0;
    }
    if (!convertToEnum(static_cast<string>(argv[3])))
    {
        return 0;
    }

    variableNames = heatSource + "-" + SurroundingT + "-" + static_cast<string>(argv[3]);
    mkdir("results", 0777);
    address = projectAddress + "/" + "results/";

    generalSetup();      // Step0
    initialConditions(); // Step1
    createFile(address); // Creating files

    for (int k = 0; k < kt; k++) // Step 10: Repeating steps 2 to 9
    {
        computeInlets(k);        // step2-part1
        liquidSurface();         // Step2-part2
        filmProperties();        // Step3
        nonDimensionalNumbers(); // Step4
        constants();             // Step6
        printData(k);            // Printing solution data
        saveData(k, address);    // Saving solution info
        temperature();           // Step7
        evaporationRate();       // Step8
        if (fluid.mLiq == 0)     // check if the liquid is ended
        {
            return 0;
        }
        fuelNewValues(); // Step9
    }

    return 0;
}

bool checkInput()
{
    if (heatSource != "constantQ" && heatSource != "variableQ")
    {
        cout << "heat source is: " << heatSource << '\n';
        cout << " heat source type is wrong!quitting...\n\n";
        return false;
    }

    if (SurroundingT != "constantT" && SurroundingT != "variableT")
    {
        cout << "surrounding T is: " << SurroundingT << '\n';
        cout << " surrounding T type is wrong!quitting...\n\n";
        return false;
    }

    return true;
}

void createFile(string &address)
{
    // creating directories
    mkdir("results/evpRate", 0777);
    mkdir("results/temperature", 0777);
    mkdir("results/nusselt", 0777);
    mkdir("results/spalding", 0777);
    mkdir("results/mass", 0777);
    mkdir("results/depth", 0777);
    mkdir("results/grashof", 0777);
    // mkdir("results/filmDensity", 0777);
    mkdir("results/rayleigh", 0777);
    mkdir("results/prandtl", 0777);
    // mkdir("results/property", 0777);
    mkdir("results/log", 0777);

    string evpRateAddress = address + "evpRate/" + variableNames;
    string temperatureAddress = address + "temperature/" + variableNames;
    string nusseltAddress = address + "nusselt/" + variableNames;
    string spaldingAddress = address + "spalding/" + variableNames;
    string massAddress = address + "mass/" + variableNames;
    string depthAddress = address + "depth/" + variableNames;
    string grashofAddress = address + "grashof/" + variableNames;
    // string filmDensityAddress = address + "filmDensity/" + variableNames;
    string rayleighAddress = address + "rayleigh/" + variableNames;
    string prandtlAddress = address + "prandtl/" + variableNames;
    // string propertyAddress = address + "property/" + variableNames;
    string logAddress = address + "log/" + variableNames;

    ofstream writeEvpRate(evpRateAddress, ios::out);         // Evaporation Rate file
    ofstream writeTemperature(temperatureAddress, ios::out); // surface Temperature file
    ofstream writeNusselt(nusseltAddress, ios::out);         // Nusselt file
    ofstream writeSpalding(spaldingAddress, ios::out);       // Spalding file
    ofstream writeMass(massAddress, ios::out);               // Liquid mass file
    ofstream writeDepth(depthAddress, ios::out);             // Liquid depth file
    ofstream writeGrashof(grashofAddress, ios::out);         // Grashof number file
    // ofstream writeFilmDensity(filmDensityAddress, ios::out); // film density file
    ofstream writePrandtl(prandtlAddress, ios::out);   // Prandtl number file
    ofstream writeRayleigh(rayleighAddress, ios::out); // Rayleigh file
    // ofstream writeProperty(propertyAddress, ios::out); // property file
    ofstream writeLog(logAddress, ios::out); // Log file

#ifdef DEBUG
    cout << " create files implemented successfully\n";
#endif
}

/* This function will save the computed data in the  file */
void saveData(int k, string &address)
{
    string evpRateAddress = address + "evpRate/" + variableNames;
    string temperatureAddress = address + "temperature/" + variableNames;
    string nusseltAddress = address + "nusselt/" + variableNames;
    string spaldingAddress = address + "spalding/" + variableNames;
    string massAddress = address + "mass/" + variableNames;
    string depthAddress = address + "depth/" + variableNames;
    string grashofAddress = address + "grashof/" + variableNames;
    // string filmDensityAddress = address + "filmDensity/" + variableNames;
    string rayleighAddress = address + "rayleigh/" + variableNames;
    string prandtlAddress = address + "prandtl/" + variableNames;
    // string propertyAddress = address + "property/" + variableNames;
    string logAddress = address + "log/" + variableNames;
    // Evaporation Rate
    ofstream writeEvpRate(evpRateAddress, ios::app);
    writeEvpRate << k * delt << "\t" << fluid.evpRate << '\n';

    // Temperature of liquid at free surface
    ofstream writeTemperature(temperatureAddress, ios::app);
    writeTemperature << k * delt << "\t" << surf.T << '\n';

    // Nusselt Number
    ofstream writeNusselt(nusseltAddress, ios::app);
    writeNusselt << k * delt << "\t" << Nu << '\n';

    // Spalding mass transfer Number
    ofstream writeSpalding(spaldingAddress, ios::app);
    writeSpalding << k * delt << "\t" << BM << '\n';

    // Grashof Number
    ofstream writeGrashof(grashofAddress, ios::app);
    writeGrashof << k * delt << "\t" << Gr << '\n';

    // Mass
    ofstream writeMass(massAddress, ios::app);
    writeMass << k * delt << "\t" << fluid.mLiq << '\n';

    // Depth
    ofstream writeDepth(depthAddress, ios::app);
    writeDepth << k * delt << "\t" << fluid.dLiq << '\n';

    // Film total density
    // ofstream writeFilmDensity(filmDensityAddress, ios::app);
    // writeFilmDensity << k * delt << "\t" << abs(film.rho - gas.rho) << '\n';

    // Rayleigh Number
    ofstream writeRayleigh(rayleighAddress, ios::app);
    writeRayleigh << k * delt << "\t" << Ra << '\n';

    // Prandtl Number
    ofstream writePrandtl(prandtlAddress, ios::app);
    writePrandtl << k * delt << "\t" << Pr << '\n';

    // Simulation Results
    ofstream writeLog(logAddress, ios::app);
    writeLog << " ********** Time step: " << k * delt << '\n'
             << " Evaporation Rate: " << fluid.evpRate << '\n'
             << " Fluid Free Surface Temperature: " << surf.T << '\n'
             << " Nusselt Number: " << Nu << '\n'
             << " Rayleigh Number: " << Ra << '\n'
             << " Grashof Number: " << Gr << '\n'
             << " Prandtl Number: " << Pr << '\n'
             << " Liquid Mass: " << fluid.mLiq << '\n'
             << " Liquid Depth: " << fluid.dLiq << '\n'
             << '\n';

#ifdef DEBUG
    cout << " save data implemented successfully\n";
#endif
}

/* STEP 0: Setting the general parameters */
void generalSetup()
{
    // Fluid Type
    // fluid.name = water; // Options: water, heptane, methanol, acetone, userFluid

    // Applying fluid and gas constant properties
    fluidConstantProperties(fluid); // DO NOT EDIT OR DELETE!!!
    gasConstantProperties(gas);     // DO NOT EDIT OR DELETE!!!

    // Heat flux

    // Computations method settings
    spaldingMethod = massFraction;        // Molar or Mass basis computation
    convectionMethod = naturalConvection; // Convection mode - Options: naturalConvection, turbulentForcedConvection

    // Temporal settings
    delt = 0.1; // second
    t = 1000.0; // second
    kt = int(t / delt) + 1;

    // Pool Characteristics
    pool.D = 1.0;                       // m
    pool.A = pi * pow(pool.D, 2) / 4.0; // m2
    pool.P = pi * pool.D;               // m
    switch (convectionMethod)
    {
    case naturalConvection:
        pool.L = pool.A / pool.P; // m
        break;
    case turbulentForcedConvection:
        pool.L = pool.D; // m
        break;
    }

#ifdef DEBUG
    cout << " general setup implemented successfully\n";
#endif
}

/* STEP 1: Initial Conditions*/
void initialConditions()
{
    // At pool
    fluid.dLiq = 0.1;                                // m
    fluid.mLiq = fluid.rhoLiq * fluid.dLiq * pool.A; // kg
    surf.T = 298.0;                                  // Kelvin
    // At surroundings
    fluid.Xsurr = 0.0; // Fluid vapor mole fraction in the surroundings
    fluid.Ysurr = 0.0; // Fluid vapor mass fraction in the surroundings

    // Gas initial conditions
    Tamb = 298.0; // Kelvin

    gas.Xsurr = 1.0 - fluid.Xsurr; // !!! DO NOT CHANGE
    gas.Ysurr = 1.0 - fluid.Ysurr; // !!! DO NOT CHANGE

#ifdef DEBUG
    cout << " initial conditions implemented successfully\n";
#endif
}

/* STEP 2: Computing the molar and mass fractions of vapor in liquid surface */
void computeInlets(int k)
{
    double transientPhase = 100.0; // time - seconds

    double qMax = 16;
    double fireTemperature = 1298; // Kelvin

    if (heatSource == "constantQ")
    {
        q = qMax;
    }
    if (heatSource == "variableQ")
    {
        if (k * delt < transientPhase)
        {
            q = k * delt * qMax / transientPhase; // kW/m2
        }
        else
        {
            q = qMax; // kW/m2
        }
    }

    if (SurroundingT == "constantT")
    {
        surr.T = fireTemperature; // Kelvin
    }
    if (SurroundingT == "variableT")
    {
        if (k * delt < transientPhase)
        {
            surr.T = 298 + k * delt * (fireTemperature - 298) / transientPhase; // Kelvin
        }
        else
        {
            surr.T = fireTemperature; // Kelvin
        }
    }

#ifdef DEBUG
    cout << " compute inlets implemented successfully\n";
#endif
}

void liquidSurface()
{
    // Vapor
    fluid.Xsurf = exp(-1.0 * (fluid.Lv * fluid.MW / R) * (1.0 / surf.T - 1.0 / fluid.Tb));         // Equation 6
    fluid.Ysurf = fluid.Xsurf / ((gas.MW / fluid.MW) + fluid.Xsurf * (1.0 - (gas.MW / fluid.MW))); // Equation 5
    // Gas
    gas.Xsurf = 1.0 - fluid.Xsurf; // mole fraction of gas
    gas.Ysurf = 1.0 - fluid.Ysurf; // mass fraction of gas

#ifdef DEBUG
    cout << " liquid surface implemented successfully\n";
#endif
}

/* STEP 3: Computing the k, mu, and Cp of film layer */
void filmProperties()
{
    film.T = surf.T + (surr.T - surf.T) / 3.0; // Equation 9 - unit: Kelvin

    fluidVariableProperties(fluid, film.T); // Computing variable properties of fluid at film temperature
    gasVariableProperties(gas, film.T);     // Computing variable properties of gas at film temperature

    fluid.Yfilm = fluid.Ysurf + (fluid.Ysurr - fluid.Ysurf) / 3.0; // Equation 10
    gas.Yfilm = 1.0 - fluid.Yfilm;                                 // Gas mass fraction at film layer

    fluid.Xfilm = fluid.Xsurf + (fluid.Xsurr - fluid.Xsurf) / 3.0; // Equation 25
    gas.Xfilm = 1.0 - fluid.Xfilm;                                 // Gas mole fraction at film layer

    film.mu = fluid.Yfilm * fluid.muVap + gas.Yfilm * gas.mu; // Equation 27

    film.k = fluid.Yfilm * fluid.kVap + gas.Yfilm * gas.k; // Equation 8

    film.MW = 1.0 / ((fluid.Xfilm / fluid.MW) + (gas.Xfilm / gas.MW)); // Equation 24

    film.rho = (Patm * film.MW) / (R * film.T); // Equation 23

    film.nu = film.mu / film.rho; // Equation 26

#ifdef DEBUG
    cout << " film properties implemented successfully\n";
#endif
}

/* STEP 4&5: Computing all non-dimensional numbers*/
void nonDimensionalNumbers()
{
    // STEP 4: Computing Grashof and Prandtl numbers
    gasVariableProperties(gas, surr.T);
    Gr = g * pow(pool.L, 3) * abs(film.rho - gas.rho) / (pow(film.nu, 2) * film.rho); // Equation 21
    Pr = gas.Cp * gas.mu / gas.k;                                                     // Equation 14
    Ra = Gr * Pr;                                                                     // Rayleigh Number

    // STEP 5: Computing Nusselt, Lewis, and Sherwood numbers
    if (film.T > surr.T) // Upper surface of heated plate or lower surface of cooled plate
    {
        if (Ra < 2e7) // Laminar flow
        {
            Nu = 0.54 * pow(Ra, 1.0 / 4.0); // Equation 18
        }
        else if (Ra > 2e7) // Turbulent flow
        {
            Nu = 0.14 * pow(Ra, 1.0 / 3.0); // Equation 19
        }
    }

    else if (film.T <= surr.T) // Upper surface of cooled plate or lower surface of heated plate
    {
        Nu = 0.27 * pow(Ra, 1.0 / 4.0); // Equation 20
    }

    switch (spaldingMethod) // Computing the Spalding mass transfer number
    {
    case massFraction:                                // In this case density of gas mixture is constant
        BM = (fluid.Ysurf - fluid.Ysurr) / gas.Ysurf; // Equation 3
        break;

    case moleFraction:                                // In this case total molar concentration and diffusivity are constant
        BM = (fluid.Xsurf - fluid.Xsurr) / gas.Xsurf; // Equation 4
        break;
    }

    Le = 1.0; // Lewis number
    Sc = 0.6; // Schmidt number
    Sh = Nu;  // Sherwood Number

#ifdef DEBUG
    cout << " non dimensional numbers implemented successfully\n";
#endif
}

/* STEP 6: This function computes the constants and Tstar and Taw */
void constants()
{
    fluidVariableProperties(fluid, surf.T);
    gasVariableProperties(gas, surf.T);

    a1 = (pool.L * fluid.eL * sigma) / (film.k * Nu); // Equation 39
    a2 = (Sh / Nu) * (fluid.Lv) / (Le * fluid.CpVap); // Equation 40

    Taw = (fluid.rhoLiq * fluid.CpLiq * pool.L * fluid.dLiq) / (film.k * Nu); // Equation 38
    switch (spaldingMethod)
    {
    case massFraction:
        ddt = (pow(fluid.Ysurf, 2) * gas.MW * fluid.Lv) / (gas.Ysurf * fluid.Xsurf * R * pow(surf.T, 2)); // Equation 42
        break;

    case moleFraction:
        ddt = (fluid.Xsurf * fluid.MW * fluid.Lv) / (gas.Xsurf * R * pow(surf.T, 2)); // Equation 43
        break;
    }

    a3 = 1.0 + (delt / Taw) * (2 * a1 * pow(surf.T, 3) + 0.5 * a2 * ddt);          // Equation 35
    a4 = 1.0 + (delt / Taw) * (a1 * pow(surf.T, 3) + 0.5 * a2 * ddt);              // Equation 36
    Tstar = a1 * pow(Tamb, 4) - a2 * log(1.0 + BM) + (q * pool.L) / (film.k * Nu); // Equation 37

#ifdef DEBUG
    cout << " constants implemented successfully\n";
#endif
}

/* STEP 7: This function computes the new temperature - temperature at new time-step */
void temperature()
{
    if (surf.T < fluid.Tb)
    {
        surfTnew = (a4 * surf.T + (delt / Taw) * Tstar) / a3; // Equation 33
        if (abs(surfTnew) > fluid.Tb)
        {
            surfTnew = fluid.Tb;
        }
    }
    else
    {
        surfTnew = fluid.Tb;
    }

#ifdef DEBUG
    cout << " temperature implemented successfully\n";
#endif
}

/* STEP 8: This function computes the evaporation rate */
void evaporationRate()
{
    if (surfTnew < fluid.Tb)
    {
        fluidVariableProperties(fluid, surf.T);
        fluid.evpRate = ((film.k * Nu) / (pool.L * fluid.CpVap * Le)) * (log(1 + BM) + 0.5 * (ddt) * (surfTnew - surf.T)); // Equation 34
    }

    else
    {
        fluid.evpRate = (q - fluid.eL * sigma * (pow(fluid.Tb, 4) - pow(Tamb, 4))) / fluid.Lv; // Equation 29
    }

    surf.T = surfTnew; // Inserting the new value

#ifdef DEBUG
    cout << " evaporation rate implemented successfully\n";
#endif
}

/* STEP 10: This function computes the new values of remaining fuel mass and new depth of fuel in the pool */
void fuelNewValues()
{
    fluid.mLiq = fluid.mLiq - delt * fluid.evpRate * pool.A; // Equation 44
    if (fluid.mLiq <= 0)
    {
        fluid.mLiq = 0.0;
        fluid.evpRate = 0.0;
    }
    fluid.dLiq = fluid.mLiq / (fluid.rhoLiq * pool.A); // Equation 45
}

/* This function will print the simulation Information in Terminal */
void printData(int k)
{
    cout << " ********** Time step: " << k * delt << '\n';
    cout << " Evaporation Rate: " << fluid.evpRate << '\n';
    cout << " Fluid Free Surface Temperature: " << surf.T << '\n';
    cout << " Nusselt Number: " << Nu << '\n';
    cout << " Rayleigh Number: " << Ra << '\n';
    cout << " Grashof Number: " << Gr << '\n';
    cout << " Prandtl Number: " << Pr << '\n';
    cout << " Liquid Mass: " << fluid.mLiq << '\n';
    cout << " Liquid Depth: " << fluid.dLiq << '\n';
    cout << '\n';

#ifdef DEBUG
    cout << " print data implemented successfully\n";
#endif
}

bool convertToEnum(string fluidName)
{
    if (fluidName == "water")
    {
        fluid.name = water;
        return true;
    }
    else if (fluidName == "heptane")
    {
        fluid.name = heptane;
        return true;
    }
    else if (fluidName == "methanol")
    {
        fluid.name = methanol;
        return true;
    }
    else if (fluidName == "acetone")
    {
        fluid.name = acetone;
        return true;
    }
    else if (fluidName == "userFluid")
    {
        fluid.name = userFluid;
        return true;
    }
    else
    {
        cout << " chosen fluid is not defined in data base! quitting...\n";
        return false;
    }
}
