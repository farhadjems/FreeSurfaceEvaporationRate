#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>

using namespace std;

/* Some Notes */
// This code computes the evaporation rate of liquid pool subjected to a constant heat fluxbased on the article 
// "Development of a numerical model for liquid pool evaporation"
// By Tarek Beji and Bart Merci
// Some notation has changed due to a week notation presented in the article
// Dalton Law is applied to compute the mixture Properties, 
// The Temperature is the same for mixture and vapor and gas in the prescribed domain 
// The Physical domain is devied into 3 sub-domians: 1.surface of liquid (surf) 2.film layer (film) 3.surrounding (surr)
// 2 substances are considered in this simulation: 1. Fluid 2. Gas
// Geometry of liquid pool is defined as a pool



/********** DEFINITION OF VARIABLES **********/

/* Definition of general constants */
double pi = 4 * atan(1); // The PI number
double sigma = 5.6704*1e-11; // Stephan-Boltzmann constant - unit: kW/(m2K4)
double R = 8.314; // Universal gas constant - unit: kj/(kmol.K)
double g = 9.81; // Gravitational acceleration - unit: m/s2  
double Patm = 101.325; // Atmospheric pressure - unit: kPa

/* Definition of non-dimensional numbers */
double Re; // Reynolds Number
double Sc; // schmidt Number
double Nu; // Nusselt Number
double Le; // Lewis Number
double BM; // Spalding mass transfer number
enum spaldingComputationMethod
{
    massFraction, moleFraction
};
spaldingComputationMethod spaldingMethod;

double Gr; // Grashof Number
double Pr; // Prandtl Number
double Ra; // Rayleigh Number
double Sh; // Sherwood Number
enum convectionComputationMethod
{
    naturalConvection, turbulentForcedConvection
};
convectionComputationMethod convectionMethod;

/* Definition of Substances */
// Fluid
struct fluidProperties
{
    double Tb; // Boiling Temperature of fluid - unit: Kelvin
    double Lv; // Latent heat of boiling - unit: kJ/kg
    double MW; // Molar weight - unit: kg/kmol
    // Liquid phase Properties
    double kLiq; // Thermal conductivity - unit: kW/(m.K)
    double CpLiq; // Thermal specific heat capacity of constant pressure - unit: kJ/(kg.K)
    double muLiq; // Viscosity - unit: Pa.s
    double rhoLiq; // Density - unit: kg/m3
    double nuLiq; // Kinematic Viscosity - unit: m3/s2
    double eL; // Emmisivity - non-dimenstional
    double alpha; // Absorptivity - non-dimenstional
    double mLiq; // Mass of liquid (fuel) - unit: kg
    double evpRate; // Evaporation rate of liquid in free surface per unit area (MLRPUA) - unit: kg/s
    double dLiq; // Liquid depth - unit: m
    // Vapor Properties
    double kVap; // Thermal conductivity - unit: kW/(m.K)
    double CpVap; // Thermal specific heat capacity of constant pressure - unit: kJ/(kg.K)
    double muVap; // Viscosity - unit: Pa.s
    double rhoVap; // Density - unit: kg/m3
    double nuVap; // Kinematic Viscosity - unit: m3/s2
    double Ysurf; // Mass fraction of vapor at surface of liquid pool
    double Xsurf; // Mole fraction of vapor at surface of liquid pool
    double Yfilm; // Mass fraction of vapor at film layer
    double Xfilm; // Mole fraction of vapor at film layer
    double Ysurr; // Mass fraction of vapor at surrounding
    double Xsurr; // Mole fraction of vapor at surrounding
};
fluidProperties fluid;
enum fluidType // The type of fluid in simulation
{
    water, heptane, acetone, methanol, userFluid,
};
fluidType fluidName; // Determines the selected liquid (fuel)

// Gas
struct gasProperties
{
    double MW; // Molar weight - unit: kg/kmol
    double k; // Thermal conductivity - unit: kW/(m.K)
    double Cp; // Thermal specific heat capacity of constant pressure - unit: kj/(kg.K)
    double mu; // Viscosity - unit: Pa.s
    double rho; // Density - unit: kg/m3
    double nu; // Kinematic Viscosity - unit: m3/s2
    double Ysurf; // Mass fraction of vapor at surface of liquid pool
    double Xsurf; // Mole fraction of vapor at surface of liquid pool
    double Yfilm; // Mass fraction of vapor at film layer
    double Xfilm; // Mole fraction of vapor at film layer
    double Ysurr; // Mass fraction of vapor at surrounding
    double Xsurr; // Mole fraction of vapor at surrounding
};
gasProperties gas;

/* Definition of domains*/
struct domainInfo
{
    double MW; // Molar weight - unit: kg/kmol
    double T; // Temperature - unit: Kelvin
    double k; // Thermal conductivity - unit: kW/(m.K)
    double Cp; // Thermal specific heat capacity of constant pressure - unit: kJ/(kg.K)
    double mu; // Viscosity - unit: Pa.s
    double rho; // Density - unit: kg/m3
    double nu; // Kinematic Viscosity - unit: m3/s2
};
domainInfo surf; // Surface of liquid
domainInfo film; // Film layer
domainInfo surr; // Surroundings

/* Definition of variable properties coefficients*/
struct vaporPropertiesInfo
{
    double a0;
    double a1;
    double a2;
    double a3;
    double a4; 
    double a5;
};
vaporPropertiesInfo k;
vaporPropertiesInfo Cp;
vaporPropertiesInfo mu;

/* Definition of Liquid Pool */
struct poolInfo
{
    double A; // Area - unit: m2
    double P; // Perimeter - unit: m
    double D; // Diameter - unit: m
    double V; // Volume - unit: m3 
    double L; // Ratio of area to perimeter - unit: m
};
poolInfo pool;

/* Definition of other variables */

// Heat flux per unit area 
double q; // Heat flux - unit: kW/m2 

// Temporal variables 
double delt; // Time-step size - unit: seconds
int kt; // Number of time steps
double t; // Total simulation time - unit: seconds

// Constants - Equations 35 - 43
double a1; // Equation 39 - unit: K-3
double a2; // Equation 40 - unit: K
double a3; // Equation 35 - unit: non-dimensional
double a4; // Equation 36 - unit: non-dimensional
double Tamb; // Ambient temperature in Radiation calculations - unit: Kelvin
double Tstar; // Equation 37 - unit: Kelvin
double Taw; // Equation 38 - unit: second
double ddt; // Equations 42 & 43 - unit: (Kelvin)-1
double surfTnew; // Temperature at new time-step - unit: Kelvin


/********** DEFINITION OF FUNCTIONS **********/
/* This function defines the constant properties of fluid*/
void fluidConstantProperties(); // Using Table 1 and Table 2 and values presented in the article 
void gasConstantProperties();
/* This function computes the variable properties of fluid using temperature */
void fluidVariablePropertiesCoefficients(double TK);
void fluidVariableProperties(double TK); // Relations used here are persented in the Appendix

void gasVariablePropertiesCoefficients(double TK);
void gasVariableProperties(double TK);

/* STEP 0: Setting the general parameters */
void generalSetup()
{
    // Fluid Type
    fluidName = water; // Options: water, heptane, methanol, acetone, userFluid
    // Applying fluid and gas constant properties
    fluidConstantProperties(); // DO NOT EDIT OR DELETE!!!
    gasConstantProperties(); // DO NOT EDIT OR DELETE!!!

    // Heat flux 
    q = 50.0; // kW/m2

    // Computations method settings    
    spaldingMethod = massFraction; // Molar or Mass basis computation
    convectionMethod = naturalConvection; // Convection mode - Options: naturalConvection, turbulentForcedConvection

    // Temporal settings
    delt = 0.1; // second
    t = 500.0; // second
    kt = int ( t/delt ) + 1;

    // Pool Characteristics
    pool.D = 0.096; // m
    pool.A = pi * pow(pool.D,2)/4.0; // m2
    pool.P = pi * pool.D; // m
    switch (convectionMethod)
    {
    case naturalConvection:
        pool.L = pool.A / pool.P; // m
        break;
    case turbulentForcedConvection:
        pool.L = pool.D; // m
        break;
    } 
} 

/* STEP 1: Initial Conditions*/
void initialConditions()
{
    // At pool
    fluid.dLiq = 0.0138; // m
    fluid.mLiq = fluid.rhoLiq * fluid.dLiq * pool.A; // kg
    surf.T = 298.0; // Kelvin
    // At surroundings
    fluid.Xsurr = 0.0; // Fluid vapor mole fraction in the surroundigns
    fluid.Ysurr = 0.0; // Fluid vapor mass fraction in the surroundigns

    // Gas initial conditions
    surr.T = 298.0; // Kelvin
    Tamb = 298.0; // Kelvin

    gas.Xsurr = 1.0 - fluid.Xsurr; // !!! DO NOT CHANGE 
    gas.Ysurr = 1.0 - fluid.Ysurr; // !!! DO NOT CHANGE 
}

/* STEP 2: Computing the molar and mass fractions of vapor in liquid surface */
void liquidSurface()
{
    // Vapor
    fluid.Xsurf = exp( -1.0 * ( fluid.Lv * fluid.MW / R ) * ( 1.0/surf.T - 1.0/fluid.Tb ) ); // Equation 6 
    fluid.Ysurf = fluid.Xsurf / ( ( gas.MW/fluid.MW ) + fluid.Xsurf * ( 1.0- ( gas.MW/fluid.MW )) ); // Equation 5 
    // Gas
    gas.Xsurf = 1.0 - fluid.Xsurf; // mole fraction of gas
    gas.Ysurf = 1.0 - fluid.Ysurf; // mass fraction of gas
}

/* STEP 3: Computing the k, mu, and Cp of film layer */
void filmProperties()
{
    film.T = surf.T + ( surr.T - surf.T )/3.0; // Equation 9 - unit: Kelvin

    fluidVariableProperties(film.T); // Computing variable properties of fluid at film temperature
    gasVariableProperties(film.T); // Computing variable properties of gas at film temperature

    fluid.Yfilm = fluid.Ysurf + ( fluid.Ysurr - fluid.Ysurf ) / 3.0; // Equation 10 
    gas.Yfilm = 1.0 - fluid.Yfilm; // Gas mass fraction at film layer
    
    fluid.Xfilm = fluid.Xsurf + ( fluid.Xsurr - fluid.Xsurf ) / 3.0; // Equation 25 
    gas.Xfilm = 1.0 - fluid.Xfilm; // Gas mole fraction at film layer
    
    film.mu = fluid.Yfilm * fluid.muVap + gas.Yfilm * gas.mu; // Equation 27

    film.k = fluid.Yfilm * fluid.kVap + gas.Yfilm * gas.k; // Equation 8
    
    film.MW = 1.0 / ( ( fluid.Xfilm/fluid.MW ) + ( gas.Xfilm/gas.MW ) ); // Equation 24

    film.rho = ( Patm * film.MW ) / ( R * film.T ); // Equation 23

    film.nu = film.mu / film.rho; // Equation 26
}

/* STEP 4&5: Computing all non-dimensional numbers*/
void nonDimensionalProperties()
{
    // STEP 4: Computing Grashof and Prandtl numbers
    gasVariableProperties(surr.T);
    Gr = g * pow(pool.L,3) * abs(film.rho - gas.rho) / ( pow(film.nu,2) * film.rho); // Equation 21 
    Pr = gas.Cp * gas.mu / gas.k; // Equation 14
    Ra = Gr * Pr; // Rayleigh Number

    // STEP 5: Computing Nusselt, Lewis, and Sherwood numbers
    if ( film.T > surr.T ) // Upper surface of heated plate or lower surface of cooled plate
    {
        if ( Ra < 2e7 ) // Laminar flow
        {
            Nu = 0.54 * pow(Ra,1.0/4.0); // Equation 18 
        }
        else if( Ra > 2e7 ) // Turbulent flow
        {
            Nu = 0.14 * pow(Ra,1.0/3.0); // Equation 19
        }
    }

    else if ( film.T <= surr.T ) // Upper surface of cooled plate or lower surface of heated plate
    {
        Nu = 0.27 * pow(Ra,1.0/4.0); // Equation 20
    }
   
    switch (spaldingMethod) // Computing the Spalding mass transfer number
    {
    case massFraction: // In this case density of gas mixture is constant
        BM = (fluid.Ysurf - fluid.Ysurr) / gas.Ysurf; // Equation 3  
        break;
    
    case moleFraction: // In this case total molar concentration and diffusivity are constant
        BM = (fluid.Xsurf - fluid.Xsurr) / gas.Xsurf; // Equation 4
        break;
    }   

    Le = 1.0; // Lewis number 
    Sc = 0.6; // Schmidt number
    Sh = Nu; // Sherwood Number
}

/* STEP 6: This function computes the constants and Tstar and Taw */
void constants()
{
    fluidVariableProperties(surf.T);
    gasVariableProperties(surf.T);

    a1 = ( pool.L * fluid.eL * sigma ) / ( film.k * Nu ); // Equation 39
    a2 = (Sh/Nu) * (fluid.Lv) / (Le * fluid.CpVap); // Equation 40
    
    Taw = (fluid.rhoLiq * fluid.CpLiq * pool.L * fluid.dLiq) / (film.k * Nu); // Equation 38
    switch (spaldingMethod)
    {
    case massFraction:
        ddt = ( pow(fluid.Ysurf,2) * gas.MW * fluid.Lv ) / ( gas.Ysurf * fluid.Xsurf * R * pow(surf.T,2) ); // Equation 42
        break;
        
    case moleFraction:
        ddt = ( fluid.Xsurf * fluid.MW * fluid.Lv ) / ( gas.Xsurf * R * pow(surf.T,2) ); // Equation 43
        break;
    }

    a3 = 1.0 + (delt/Taw) * ( 2 * a1 * pow(surf.T,3) + 0.5 * a2 * ddt  ); // Equation 35
    a4 = 1.0 + (delt/Taw) * ( a1 * pow(surf.T,3) + 0.5 * a2 * ddt  ); // Equation 36
    Tstar = a1 * pow(Tamb,4) - a2 * log(1.0+BM) + ( q * pool.L ) / ( film.k * Nu ); // Equation 37
}

/* STEP 7: This function computes the new temperature - temperatue at new time-step */
void temperature()
{
    if (surf.T < fluid.Tb)
    {
        surfTnew = ( a4 * surf.T + (delt/Taw) * Tstar ) / a3; // Equation 33
        if ( abs(surfTnew) > fluid.Tb )
        {
            surfTnew = fluid.Tb;
        }
    }
    else
    {
        surfTnew = fluid.Tb;
    }
}

/* STEP 8: This function computes the evaporation rate */
void evaporationRate ()
{
    if ( surfTnew < fluid.Tb )
    {
        fluidVariableProperties(surf.T);
        fluid.evpRate = ( ( film.k * Nu ) / ( pool.L * fluid.CpVap * Le ) ) * ( log(1+BM) + 0.5 * (ddt) *(surfTnew - surf.T) ); // Equation 34
    }

    else
    {
        fluid.evpRate = ( q - fluid.eL * sigma * (pow(fluid.Tb,4) - pow(Tamb,4)) ) / fluid.Lv; // Equation 29
    }
    
    surf.T = surfTnew; // Inserting the new value
}

/* STEP 10: This function computes the new values of remaining fuel mass and new depth of fuel in the pool */
void fuelNewValues()
{
    fluid.mLiq = fluid.mLiq - delt * fluid.evpRate * pool.A; // Equation 44
    if (fluid.mLiq <= 0)
    {
        fluid.mLiq = 0.0;
        fluid.evpRate = 0.0;
        surf.T = 298.0;
    }
    fluid.dLiq = fluid.mLiq / (fluid.rhoLiq * pool.A); // Equation 45
}

/* This function will print the simulation Information in Terminal */
void printData(int k)
{
    cout << " ********** Time step: " << k * delt << endl ;
    cout << " Evaporation Rate: " << fluid.evpRate << endl;
    cout << " Fluid Free Surface Temperature: " << surf.T << endl;
    cout << " Nusselt Number: " << Nu << endl;
    cout << " Rayleigh Number: " << Ra << endl;
    cout << " Grashof Number: " << Gr << endl;
    cout << " Prandtl Number: " << Pr << endl;
    cout << " Liquid Mass: " << fluid.mLiq << endl; 
    cout << " Liquid Depth: " << fluid.dLiq << endl ;
    cout << endl; 
}

/* Creating files for saving th simulation data */

/* This Finction Creates the output files*/
void createFile()
{
    string evpRateAddress = "/home/farhad/Codes/Cpp/evpRate/version1/results/evpRate/"+ to_string(int(q)) + "kw-" + to_string(int(surr.T)) + "K";
    string temperatureAddress = "/home/farhad/Codes/Cpp/evpRate/version1/results/temperature/" + to_string(int(q))+ "kw-" + to_string(int(surr.T)) + "K"; 
    string nusseltAddress = "/home/farhad/Codes/Cpp/evpRate/version1/results/nusselt/"+ to_string(int(q))+ "kw-" + to_string(int(surr.T)) + "K";
    string spaldingAddress = "/home/farhad/Codes/Cpp/evpRate/results/spalding/"+ to_string(int(q))+  "kw-" +to_string(int(surr.T)) + "K";
    string massAddress = "/home/farhad/Codes/Cpp/evpRate/results/mass/"+ to_string(int(q))+  "kw-" + to_string(int(surr.T)) + "K";
    string depthAddress = "/home/farhad/Codes/Cpp/evpRate/results/depth/"+ to_string(int(q))+  "kw-" + to_string(int(surr.T)) + "K";
    string grashofAddress = "/home/farhad/Codes/Cpp/evpRate/results/grashof/"+ to_string(int(q))+  "kw-" + to_string(int(surr.T)) + "K";
    string filmDensityAddress = "/home/farhad/Codes/Cpp/evpRate/results/filmDensity/"+ to_string(int(q))+  "kw-" + to_string(int(surr.T)) + "K";  
    // string rayleighAddress = "/home/farhad/Codes/Cpp/evpRate/results/Rayleigh";
    // string prandtlAddress = "/home/farhad/Codes/Cpp/evpRate/results/Prandtl";
    // string propertyAddress = "/home/farhad/Codes/Cpp/evpRate/results/property";
    // string logAddress = "/home/farhad/Codes/Cpp/evpRate/results/Log";
    
    ofstream writeEvpRate(evpRateAddress, ios::out); // Evaporation Rate file
    ofstream writeTemperature(temperatureAddress, ios::out); // surface Temperature file
    ofstream writeNusselt(nusseltAddress, ios::out); // Nusselt file
    ofstream writeSpalding(spaldingAddress, ios::out); // Spalding file
    ofstream writeMass(massAddress, ios::out); // Liquid mass file
    ofstream writeDepth(depthAddress, ios::out); // Liquid depth file
    ofstream writeGrashof(grashofAddress, ios::out); // Grashof number file
    ofstream writeFilmDensity(filmDensityAddress, ios::out); // film density file
    
    /*
        ofstream writePrandtl(prandtlAddress, ios::out); // Prandtl number file
        ofstream writeRayleigh(rayleighAddress, ios::out); // Rayleigh file
        ofstream writeProperty(propertyAddress, ios::out); // property file    
        ofstream writeLog(logAddress, ios::out); // Log file 
    */
}

/* This function will save the computed data in the  file */
void saveData(int k)
{
    string evpRateAddress = "/home/farhad/Codes/Cpp/evpRate/version1/results/evpRate/"+ to_string(int(q)) + "kw-" + to_string(int(surr.T)) + "K";
    string temperatureAddress = "/home/farhad/Codes/Cpp/evpRate/version1/results/temperature/" + to_string(int(q))+ "kw-" + to_string(int(surr.T)) + "K"; 
    string nusseltAddress = "/home/farhad/Codes/Cpp/evpRate/version1/results/nusselt/"+ to_string(int(q))+ "kw-" + to_string(int(surr.T)) + "K";
    string spaldingAddress = "/home/farhad/Codes/Cpp/evpRate/results/spalding/"+ to_string(int(q))+  "kw-" +to_string(int(surr.T)) + "K";
    string massAddress = "/home/farhad/Codes/Cpp/evpRate/results/mass/"+ to_string(int(q))+  "kw-" + to_string(int(surr.T)) + "K";
    string depthAddress = "/home/farhad/Codes/Cpp/evpRate/results/depth/"+ to_string(int(q))+  "kw-" + to_string(int(surr.T)) + "K";
    string grashofAddress = "/home/farhad/Codes/Cpp/evpRate/results/grashof/"+ to_string(int(q))+  "kw-" + to_string(int(surr.T)) + "K";
    string filmDensityAddress = "/home/farhad/Codes/Cpp/evpRate/results/filmDensity/"+ to_string(int(q))+  "kw-" + to_string(int(surr.T)) + "K"; 
    // string rayleighAddress = "/home/farhad/Codes/Cpp/evpRate/results/Rayleigh";
    // string prandtlAddress = "/home/farhad/Codes/Cpp/evpRate/results/Prandtl";
    // string propertyAddress = "/home/farhad/Codes/Cpp/evpRate/results/property";
    // string logAddress = "/home/farhad/Codes/Cpp/evpRate/results/Log";
    // Evaporation Rate
    ofstream writeEvpRate(evpRateAddress, ios::app);
    writeEvpRate << k * delt << "\t" << fluid.evpRate << endl; 
    
    // Temperature of liquid at free surface
    ofstream writeTemperature(temperatureAddress, ios::app);
    writeTemperature << k * delt << "\t" << surf.T << endl;     
    
    // Nusselt Number
    ofstream writeNusselt(nusseltAddress, ios::app);
    writeNusselt << k * delt << "\t" << Nu << endl;    
    
    // Spalding mass tansfer Number
    ofstream writeSpalding(spaldingAddress, ios::app);
    writeSpalding << k * delt << "\t" << BM << endl; 

    // Grashof Number
    ofstream writeGrashof(grashofAddress, ios::app);
    writeGrashof << k * delt << "\t" << Gr << endl; 

    //Mass
    ofstream writeMass(massAddress, ios::app);
    writeMass << k * delt << "\t" << fluid.mLiq << endl;   
    
    //Depth
    ofstream writeDepth(depthAddress, ios::app);
    writeDepth << k * delt << "\t" << fluid.dLiq << endl;   
    
    // Film total density 
    ofstream writeFilmDensity(filmDensityAddress, ios::app);
    writeFilmDensity << k * delt << "\t" << abs(film.rho - gas.rho) << endl; 

    // // Rayleigh Number
    // ofstream writeRayleigh(rayleighAddress, ios::app);
    // writeRayleigh << k * delt << "\t" << Ra << endl;   
    
    // // Prandtl Number
    // ofstream writePrandtl(prandtlAddress, ios::app);
    // writePrandtl << k * delt << "\t" << Pr << endl; 

    // // fluid property
    // gasVariableProperties(surf.T);
    // ofstream writeProperty(propertyAddress, ios::app);
    // writeProperty << surf.T << "\t" << gas.k << endl; 

    // // Simulation Results
    // ofstream writeLog(logAddress, ios::app);
    // writeLog << " ********** Time step: " << k * delt << endl ;
    // writeLog << " Evaporation Rate: " << fluid.evpRate << endl;
    // writeLog << " Fluid Free Surface Temperature: " << surf.T << endl;
    // writeLog << " Nusselt Number: " << Nu << endl;
    // writeLog << " Rayleigh Number: " << Ra <<endl;
    // writeLog << " Grashof Number: " << Gr << endl;
    // writeLog << " Prandtl Number: " << Pr << endl;
    // writeLog << " Liquid Mass: " << fluid.mLiq << endl; 
    // writeLog << " Liquid Depth: " << fluid.dLiq << endl ;
    // writeLog << endl;
}

/* The below function is the code algorithm. this algorithm is based on the presented algorithm in page 52 */
void algorithm()
{
    generalSetup(); // Step0
    initialConditions(); // Step1
    createFile();
    for (int k = 0; k < kt; k++) // Step 10: Repeating steps 2 to 9
    {
        liquidSurface(); // Step2
        filmProperties(); // Step3
        nonDimensionalProperties(); // Step4
        constants(); // Step6
        printData(k); // Printing solution data 
        saveData(k); // Saving solution info
        temperature(); // Step7
        evaporationRate(); // Step8
        if (fluid.mLiq == 0)
        {
            return;
        }
        
        fuelNewValues(); // Step9     
    }
}

/* Main function of code. DO NOT CHANGE IT!!! */
int main()
{
    algorithm(); // This function contains the alorithm of code

    return 0;
}


/* Definition of Property functions*/
// Fluid Properties
void fluidConstantProperties()
{
    // Table 1 - Liquid Properties 
    switch (fluidName)
    {
    case water: //H2O
        fluid.MW = 18.0; // unit: kg/kmol
        fluid.rhoLiq = 1000.0; // unit: kg/m3
        fluid.CpLiq = 4.18; // unit: kj/(kg.k)
        fluid.Lv = 2260.0; // unit: kj/kg
        fluid.eL = 0.95; // non-dimensional
        fluid.Tb = 373.0; // unit: Kelvin  
        fluid.alpha = 1.0; // non-dimensional
        break;

    case heptane: // C7H16
        fluid.MW = 100.0; // unit: kg/kmol
        fluid.rhoLiq = 675.0; // unit: kg/m3
        fluid.CpLiq = 2.24; // unit: kj/(kg.k)
        fluid.Lv = 317.0; // unit: kj/kg
        fluid.eL = 0.95; // non-dimensional
        fluid.Tb = 371.5; // unit: Kelvin
        fluid.alpha = 1.0; // non-dimensional
        break;

    case acetone: // C3H6O
        fluid.MW = 58.0; // unit: kg/kmol
        fluid.rhoLiq = 791.0; // unit: kg/m3
        fluid.CpLiq = 2.13; // unit: kj/(kg.k)
        fluid.Lv = 501.0; // unit: kj/kg
        fluid.eL = 0.95; // non-dimensional
        fluid.Tb = 329.3; // unit: Kelvin
        fluid.alpha = 1.0; // non-dimensional
        break;

    case methanol: // CH3OH
        fluid.MW = 32.0; // unit: kg/kmol
        fluid.rhoLiq = 796.0; // unit: kg/m3
        fluid.CpLiq = 2.48; // unit: kj/(kg.k)
        fluid.Lv = 1099.0; // unit: kj/kg
        fluid.eL = 0.95; // non-dimensional
        fluid.Tb = 337.8; // unit: Kelvin
        fluid.alpha = 1.0; // non-dimensional
        break;

    case userFluid: // user-defined material. enter the properties and change the variable "liquidName" value to "user"
        fluid.MW = 1.0; // unit: kg/kmol
        fluid.rhoLiq = 1.0; // unit: kg/m3
        fluid.CpLiq = 1.0; // unit: kj/(kg.k)
        fluid.Lv = 1.0; // unit: kj/kg
        fluid.eL = 1.0; // non-dimensional
        fluid.Tb = 273.0; // unit: Kelvin
        fluid.alpha = 1.0; // non-dimensional
        break;
    }
}

void fluidVariablePropertiesCoefficients(double TK)
{
    double TC = TK - 273.0;
    switch (fluidName)
    {
    case water:
        // vap.Cp coefficients
        Cp.a0 = 2.67703787e1;
        Cp.a1 = 2.97318329e-3;
        Cp.a2 = -7.73769690e-7;
        Cp.a3 = 9.44336689e-11;
        Cp.a4 = -4.26900959e-15;
        Cp.a5 = 0;
        // vap.k and vap.mu coefficients
        if ( TK >= 295 && TK <= 1000 )
        {
            k.a0 = 0.15541443e1;
            k.a1 = 0.66106305e2;
            k.a2 = 0.55968860e4;
            k.a3 = -0.39259598e1;
            k.a4 = 0.0;
            k.a5 = 0.0;

            mu.a0 = 0.78387780e0;
            mu.a1 = -0.38260408e3;
            mu.a2 = 0.49040158e5;
            mu.a3 = 0.85222785e0;
            mu.a4 = 0.0;
            mu.a5 = 0.0;
        }
        if ( TK > 1000 && TK <= 5000 )
        {
            k.a0 = 0.79349503e0;
            k.a1 = -0.13340063e4;
            k.a2 = 0.378664327e6;
            k.a3 = 0.23591474e1;
            k.a4 = 0.0;
            k.a5 = 0.0;

            mu.a0 = 0.50714993e0;
            mu.a1 = -0.68966913e3;
            mu.a2 = 0.874547505e5;
            mu.a3 = 0.30285155e1;
            mu.a4 = 0.0;
            mu.a5 = 0.0;
        }
        break;
    
    case heptane:
        // vap.Cp coefficients
        Cp.a0 = -1.3934e-1;
        Cp.a1 = 2.6295e-3;
        Cp.a2 = 1.15e-5;
        Cp.a3 = -1.7062e-7;
        Cp.a4 = -1.205e-9;
        Cp.a5 = 1.5923e-11;
        // vap.k coefficients
        k.a0 = -4.6160e0;
        k.a1 = 7.7908e-3;
        k.a2 = -1.5605e-5;
        k.a3 = 1.158e-7;
        k.a4 = -1.4660e-9;
        k.a5 = 5.4467e-12;
        // vap.mu coefficients
        mu.a0 = 3.8877;
        mu.a1 = 6.1806e-3;
        mu.a2 = -3.0481e-5;
        mu.a3 = 9.7703e-8;
        mu.a4 = 0.0;
        mu.a5 = 0.0;
        break;

    case acetone:
        // vap.Cp coefficients
        Cp.a0 = 1.9452e-1;
        Cp.a1 = 2.2900e-3;
        Cp.a2 = -8.6330e-7;
        Cp.a3 = -2.0672e-8;
        Cp.a4 = 1.9250e-10;
        Cp.a5 = -6.5969e-13;
        // vap.k coefficients
        k.a0 = -4.6467;
        k.a1 = 7.1871e-3;
        k.a2 = -2.0976e-5;
        k.a3 = 4.6070e-7;
        k.a4 = -5.4286e-9;
        k.a5 = 1.9213e-11;
        // vap.mu coefficients
        mu.a0 = 4.3519e0;
        mu.a1 = 2.6016e-3;
        mu.a2 = -3.5279e-6;
        mu.a3 = 2.1574e-7;
        mu.a4 = -3.5172e-9;
        mu.a5 = 1.3856e-11;
        break;
        
    case methanol:
        // vap.Cp coefficients
        Cp.a0 = 0.31404e0;
        Cp.a1 = 5.01e-3;
        Cp.a2 = -1.7365e-4;
        Cp.a3 = 3.4375e-6;
        Cp.a4 = -2.5001e-8;
        Cp.a5 = 5.9705e-11;
        // vap.k coefficients
        k.a0 = -6.7646e0;
        k.a1 = 8.9038e-3;
        k.a2 = -3.5238e-5;
        k.a3 = 1.4232e-7;
        k.a4 = -4.003e-10;
        k.a5 = 5.0729e-13;
        // vap.mu coefficients
        mu.a0 = 4.4794e0;
        mu.a1 = 3.3930e-3;
        mu.a2 = 4.0163e-6;
        mu.a3 = -7.3376e-8;
        mu.a4 = 2.4641e-10;
        mu.a5 = -2.1042e-13;
        break; 

    case userFluid: // Enter tthe constants here
        // vap.Cp coefficients
        Cp.a0 = 1.0;
        Cp.a1 = 1.0;
        Cp.a2 = 1.0;
        Cp.a3 = 1.0;
        Cp.a4 = 1.0;
        Cp.a5 = 1.0;
        // vap.k coefficients
        k.a0 = 1.0;
        k.a1 = 1.0;
        k.a2 = 1.0;
        k.a3 = 1.0;
        k.a4 = 1.0;
        k.a5 = 1.0;
        // vap.mu coefficients
        mu.a0 = 1.0;
        mu.a1 = 1.0;
        mu.a2 = 1.0;
        mu.a3 = 1.0;
        mu.a4 = 1.0;
        mu.a5 = 1.0;
        break;
    }  
}

void fluidVariableProperties(double TK)
{
    double TC = TK - 273.0;
    fluidVariablePropertiesCoefficients(TK);
    switch (fluidName)
    {
    case water:
        fluid.CpVap = ( Cp.a0 + Cp.a1 * TK + Cp.a2 * pow(TK,2) + Cp.a3 * pow(TK,3) + Cp.a4 * pow(TK,4) ) * R / fluid.MW; // Equation A1
        fluid.kVap = exp( k.a0 * log(TK) + k.a1 / TK + k.a2 / pow(TK,2) + k.a3 ) * 1e-7; // Equation A2
        fluid.muVap = exp( mu.a0 * log(TK) + mu.a1 / TK + mu.a2 / pow(TK,2) + mu.a3 ) * 1e-7; // Equation A3
        break;
    
    case (heptane):
        fluid.CpVap = exp( Cp.a0 + Cp.a1 * TC + Cp.a2 * pow(TC,2) + Cp.a3 * pow(TC,3) + Cp.a4 * pow(TC,4) + Cp.a5 * pow(TC,5) ); // Equation A3
        fluid.kVap = exp( k.a0 + k.a1 * TC + k.a2 * pow(TC,2) + k.a3 * pow(TC,3) + k.a4 * pow(TC,4) + k.a5 * pow(TC,5) )* 1e-3; // Equation A3
        fluid.muVap = exp( mu.a0 + mu.a1 * TC + mu.a2 * pow(TC,2) + mu.a3 * pow(TC,3) + mu.a4 * pow(TC,4) + mu.a5 * pow(TC,5) ) * 1e-7; // Equation A3
        break;

    case (acetone):
        fluid.CpVap = exp( Cp.a0 + Cp.a1 * TC + Cp.a2 * pow(TC,2) + Cp.a3 * pow(TC,3) + Cp.a4 * pow(TC,4) + Cp.a5 * pow(TC,5) ); // Equation A3
        fluid.kVap = exp( k.a0 + k.a1 * TC + k.a2 * pow(TC,2) + k.a3 * pow(TC,3) + k.a4 * pow(TC,4) + k.a5 * pow(TC,5) )* 1e-3; // Equation A3
        fluid.muVap = exp( mu.a0 + mu.a1 * TC + mu.a2 * pow(TC,2) + mu.a3 * pow(TC,3) + mu.a4 * pow(TC,4) + mu.a5 * pow(TC,5) ) * 1e-7; // Equation A3
        break;

    case (methanol):        
        fluid.CpVap = exp( Cp.a0 + Cp.a1 * TC + Cp.a2 * pow(TC,2) + Cp.a3 * pow(TC,3) + Cp.a4 * pow(TC,4) + Cp.a5 * pow(TC,5) ); // Equation A3
        fluid.kVap = exp( k.a0 + k.a1 * TC + k.a2 * pow(TC,2) + k.a3 * pow(TC,3) + k.a4 * pow(TC,4) + k.a5 * pow(TC,5) )* 1e-3; // Equation A3
        fluid.muVap = exp( mu.a0 + mu.a1 * TC + mu.a2 * pow(TC,2) + mu.a3 * pow(TC,3) + mu.a4 * pow(TC,4) + mu.a5 * pow(TC,5) ) * 1e-7; // Equation A3
        break;

    case userFluid:
        fluid.CpVap = 1.0;
        fluid.kVap = 1.0;
        fluid.muVap = 1.0;
        break;
    }
}

// Air Properties
void gasConstantProperties()
{
    gas.MW = 29.0; // unit: kg/kmol

}

void gasVariablePropertiesCoefficients(double TK)
{
    double TC = TK - 273.0;
    // Cp
    Cp.a0 = 7.5512e-3;
    Cp.a1 = 9.8834e-6;
    Cp.a2 = 5.0936e-7;
    Cp.a3 = -5.3667e-10;
    Cp.a4 = 2.1073e-13;
    Cp.a5 = -2.3807e-17;
    // k
    k.a0 = 3.1349e0;
    k.a1 = 3.5533e-3;
    k.a2 = -5.2910e-6;
    k.a3 = 4.4429e-9;
    k.a4 = -1.7077e-12;
    k.a5 = 2.4576e-16;
    // mu
    mu.a0 = 5.0907e0;
    mu.a1 = 3.2040e-3;
    mu.a2 = -4.7190e-6;
    mu.a3 = 3.8481e-9;
    mu.a4 = -1.4484e-12;
    mu.a5 = 2.0023e-16;
}

void gasVariableProperties(double TK)
{
    double TC = TK - 273.0; // Input Temperature - Celsius
    gasVariablePropertiesCoefficients(TK);
    gas.Cp = exp( Cp.a0 + Cp.a1 * TC + Cp.a2 * pow(TC,2) + Cp.a3 * pow(TC,3) + Cp.a4 * pow(TC,4) + Cp.a5 * pow(TC,5) ); // Equation A3
    gas.k = exp( k.a0 + k.a1 * TC + k.a2 * pow(TC,2) + k.a3 * pow(TC,3) + k.a4 * pow(TC,4) + k.a5 * pow(TC,5) ) * 1e-6; // Equation A3
    gas.mu = exp( mu.a0 + mu.a1 * TC + mu.a2 * pow(TC,2) + mu.a3 * pow(TC,3) + mu.a4 * pow(TC,4) + mu.a5 * pow(TC,5) ) * 1e-7; // Equation A3
    gas.rho = ( Patm * gas.MW ) / ( R * TK ); 
}

