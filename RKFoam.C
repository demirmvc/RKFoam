#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "OFstream.H"
#include "RKHelperFunctions.H"  // Include before main

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readRKSchemeSettings.H"
    #include "createFvOptions.H"

    turbulence->validate();

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

        // Initialize U0 and phi0 at the beginning of each time step
        volVectorField U0(U.name(), U);

        surfaceScalarField phi0(phi.name(), phi);

        // Runge-Kutta stages loop
        for (label stage = 0; stage < nStages; stage++)
        {
            #include "solveRKStage.H"
        }

        // Update fields after all stages
        #include "updateFields.H"

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
