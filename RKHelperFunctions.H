#ifndef RKHelperFunctions_H
#define RKHelperFunctions_H

#include "fvCFD.H"

// Function declaration
surfaceScalarField computeFlux
(
    const volVectorField& HbyA,
    const volScalarField& rAU,
    const volVectorField& U,
    const volScalarField& p,
    const surfaceScalarField& phi,
    const word& method
);

// Function definition
inline surfaceScalarField computeFlux
(
    const volVectorField& HbyA,
    const volScalarField& rAU,
    const volVectorField& U,
    const volScalarField& p,
    const surfaceScalarField& phi,
    const word& method
)
{
    const fvMesh& mesh = U.mesh();

    if (method == "Rhie-Chow")
    {
        return fvc::flux(HbyA)
             + (fvc::interpolate(rAU * (fvc::grad(p) - fvc::grad(p, "cell"))) & mesh.Sf());
    }
    else if (method == "OpenFOAM")
    {
        surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
        surfaceScalarField phiCorr
        (
            "phiCorr",
            fvc::interpolate(rAU) * fvc::ddtCorr(U, phi)
        );
        return phiHbyA + phiCorr;
    }
    else if (method == "DAlessandro")
    {
        return fvc::flux(HbyA)
             + fvc::interpolate(rAU) * fvc::ddtCorr(U, phi)
             - fvc::snGrad(p) * fvc::interpolate(rAU) * mesh.magSf();
    }
    else
    {
        FatalErrorIn("computeFlux")
            << "Unknown velocity interpolation method: " << method << endl
            << abort(FatalError);
        return surfaceScalarField::null();
    }
}

#endif
