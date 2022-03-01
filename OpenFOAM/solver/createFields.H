#include "createRDeltaT.H"
IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
IOdictionary controlDict
    (
        IOobject
        (
            "controlDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );  
//scalar startTime2 = readScalar(controlDict.subDict("kcTimeControl").lookup("startTime2"));     
//scalar endTime2 = readScalar(controlDict.subDict("kcTimeControl").lookup("endTime2"));     


//read Turbulent Schmidt number    
dimensionedScalar Sct
    (
        transportProperties.lookup("Sct")
    );   
    
//read mass diffusivity of O2   
dimensionedScalar DAB
    (
        transportProperties.lookup("DAB")
    );   
    
Info<< "Reading field CO2\n" << endl;
//
volScalarField CO2
(
    IOobject
    (
        "CO2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField DABt
(
    IOobject
    (
        "DABt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("DABt", dimensionSet (0,2,-1,0,0), 0)
);
volScalarField DABEff
(
    IOobject
    (
        "DABEff",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    DABt + DAB
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Molar flux of O2
//heat flux 
    volScalarField NO2
    (
        IOobject
        (
            "NO2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
         mesh,
         dimensionedScalar("NO2", dimensionSet (0,-2,-1,0,1), 0)
    );

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, laminarTransport)
);

//instantaneous area-averaged mass-transfer coefficient
    uniformDimensionedScalarField kcInst
    (
        IOobject
        (
            "kcInst",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar ("kcInst",dimLength/dimTime,0)
    );

//time_average mass-transfer coefficient
    uniformDimensionedScalarField kc
    (
        IOobject
        (
            "kc",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar ("kc",dimLength/dimTime,0)
    );

    uniformDimensionedScalarField CO2b
    (
        IOobject
        (
            "CO2b",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar ("CO2b",dimMoles/dimLength/dimLength/dimLength,0)
    );
    
// Fields values of previous variables: to be used with probeLocations
//instantaneous area-averaged mass-transfer coefficient
    volScalarField kcInstF
    (
        IOobject
        (
            "kcInstF",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar ("kcInstF",dimLength/dimTime,0)
    );

//time_average mass-transfer coefficient
    volScalarField kcF
    (
        IOobject
        (
            "kcF",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,     
        dimensionedScalar ("kcF",dimLength/dimTime,0)
    );

    volScalarField CO2bF
    (
        IOobject
        (
            "CO2bF",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,     
        dimensionedScalar ("CO2bF",dimMoles/dimLength/dimLength/dimLength,0)
    );

#include "createMRF.H"
#include "createFvOptions.H"
