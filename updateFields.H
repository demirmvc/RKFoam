// Final update after all stages
U = U0;
phi = phi0;

for (label stage = 0; stage < nStages; stage++)
{
    U.internalField() += runTime.deltaTValue() * b[stage];
    phi.internalField() += runTime.deltaTValue() * b[stage];
}

U.correctBoundaryConditions();
