// proxy for Convex2
// Contains all sub-routines, but different main-file
// 
// Intended for quick flux-calculations

#include "../ThermalInertiaOnlyConvex.h"
#include<iostream>
using namespace std;

int main()
{
    char shapefilename[] = "s:/spheres/cube4.obj.convex";
    ConvexFile shape (shapefilename);

    const double H = 17.365; // pV = 0.2 --> D=1km
    const double G = 0.15;
    const double emissivity = 0.9;
    const double pV = 0.2;
    const double ecl_lambda = 0;
    const double ecl_beta   = 89.999;
    const double period_h = 3;
    const double JD0 = 0;
    SpinState axis(ecl_lambda, ecl_beta, period_h, JD0);

    TriangulatedConvex aster(shape, axis, H, G, pV, emissivity);

    const double ThermalParameter = 0;
    const int nTime = 300;
    const int nZ = 25;
    const double zMax = 6;
    const double accuracyGoal = 1e-3;
    ThermalInertiaOnlyConvex model(aster, ThermalParameter, nTime, nZ, zMax, accuracyGoal);

    const vector3 A2Sun(1.,1e-4,0);
    //const vector3 A2Earth(1e-4,1.,0);
    const vector3 A2Earth(1,0,0);
    const double lambdaMu = 10;
    const double rAU = 1;
    const double deltaKM = AsteroidConstants::AU;
    double flux = aster.ThermalFlux(A2Sun, A2Earth, lambdaMu, rAU, deltaKM, model);
    cout<<flux<<endl;


    return 0;
};