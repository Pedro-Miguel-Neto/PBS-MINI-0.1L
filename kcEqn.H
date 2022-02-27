Info<< "Compute molar flux of O2 in the patches" << endl;

//O2 Molar flux
surfaceScalarField sNO2 = fvc::interpolate(DABEff)*fvc::snGrad(CO2);

//Local molar flux on the patch "top"
label patchTop = mesh.boundaryMesh().findPatchID("top");

//
volScalarField::Boundary& NO2Bf = NO2.boundaryFieldRef();

const surfaceScalarField::Boundary& sNO2Bf = sNO2.boundaryField();

NO2Bf[patchTop] = sNO2Bf[patchTop];
//


dimensionedScalar patchArea ("patchArea", dimArea, 0.0);
dimensionedScalar AvTopNO2 ("AvTopNO2", dimMoles/dimArea/dimTime, 0.0);
dimensionedScalar AvTopCO2 ("AvTopCO2", dimMoles/dimVolume, 0.0);

patchArea.value()  = gSum(mesh.magSf().boundaryField()[patchTop]);

AvTopNO2.value() = gSum
(
    mesh.magSf().boundaryField()[patchTop]*
    NO2.boundaryField()[patchTop]
)/patchArea.value();

AvTopCO2.value() = gSum
(
    mesh.magSf().boundaryField()[patchTop]*
    CO2.boundaryField()[patchTop]
)/patchArea.value();

CO2b = CO2.weightedAverage(mesh.V());
//instantaneous area-averaged mass-transfer coefficient
kcInst = AvTopNO2/(AvTopCO2-CO2b); 

scalar startTime2 = readScalar(controlDict.subDict("kcTimeControl").lookup("startTime2"));     
scalar endTime2 = readScalar(controlDict.subDict("kcTimeControl").lookup("endTime2"));  

//time_average mass-transfer coefficient
if (runTime.value() < startTime2) {kc.value() = 0;}
if (runTime.value() > startTime2 && runTime.value() < endTime2)
{
  kc.value() = 
    (
      kc.value()*(runTime.value()-startTime2)  
     +kcInst.value()*(runTime.deltaT().value())
    )
   /
   (runTime.value()-startTime2+runTime.deltaT().value());
}
Info<< "CO2.weightedAverage = " << CO2.weightedAverage(mesh.V()) <<
endl;
Info<< "AvTopNO2 = " << AvTopNO2 << endl;
Info<< "kc = " << kc << endl;
Info<< "kcInst = " << kcInst << endl;

CO2bF = CO2b;
kcF = kc;
kcInstF = kcInst;
