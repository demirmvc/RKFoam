// Store previous stage values
volVectorField U_prev = U;
surfaceScalarField phi_prev = phi;

// Outer iteration loop
for (label outer = 0; outer < nOuterIterations; outer++)
{
    // Momentum predictor
    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));
        fvOptions.correct(U);
    }

    // PISO loop
    for (label corr = 0; corr < nCorr; corr++)
    {
        volScalarField rAU(1.0/UEqn.A());
        volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            computeFlux(HbyA, rAU, U, p, velocityInterpolationMethod)
        );

        // Non-orthogonal pressure corrector loop
        for (label nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
            );

            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve();

            if (nonOrth == nNonOrthCorr)
            {
                phi = phiHbyA - pEqn.flux();
            }
        }

        #include "continuityErrs.H"

        // Explicitly relax pressure for momentum corrector
        p.relax();

        // Momentum corrector
        U = HbyA - rAU*fvc::grad(p);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}

// Update RK stage values
U = U0;
phi = phi0;

for (label stage = 0; stage < nStages; stage++)
{
U.internalField() += runTime.deltaTValue() * A[stage][i] * deltaUList[i].internalField();
phi.internalField() += runTime.deltaTValue() * A[stage][i] * deltaPhiList[i].internalField();
}
