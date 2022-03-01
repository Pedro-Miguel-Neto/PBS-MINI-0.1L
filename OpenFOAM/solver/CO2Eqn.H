// Solve the conservation equation of [O2]
if (turbulence.valid() ) 
  {
    DABt = turbulence->nut()/Sct;
    DABt.correctBoundaryConditions();
  }
else
  {
    DABt = DABt*0;  
  }

DABEff == DABt + DAB;

fvScalarMatrix CO2Eqn
(
    fvm::ddt(CO2)
  + fvm::div(phi, CO2) ==
  fvm::laplacian(DABEff,CO2) 
);

CO2Eqn.relax();
CO2Eqn.solve();



