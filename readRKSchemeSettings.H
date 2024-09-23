Info << "Reading RK scheme settings\n" << endl;

IOdictionary rkDict
(
    IOobject
    (
        "RKProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const word rkScheme(rkDict.lookup("RKScheme"));
const label nStages(readLabel(rkDict.lookup("nStages")));
const label nOuterIterations(readLabel(rkDict.lookup("nOuterIterations")));
const label nCorr(readLabel(rkDict.lookup("nCorr")));
const label nNonOrthCorr(readLabel(rkDict.lookup("nNonOrthCorr")));
const word velocityInterpolationMethod(rkDict.lookup("velocityInterpolationMethod"));

// Read Butcher tableau
List<List<scalar>> A(nStages);

for (label i = 0; i < nStages; i++)
{
    List<scalar> Ai;
    rkDict.lookup("A" + Foam::name(i + 1)) >> Ai;
    A[i] = Ai;
}

scalarField b(nStages);
rkDict.lookup("b") >> b;

scalarField c(nStages);
rkDict.lookup("c") >> c;
