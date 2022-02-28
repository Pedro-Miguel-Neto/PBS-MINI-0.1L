# PBS-MINI-0.1L
A research on the hydrodynamics and aeration of the PBS-Biotech MINI 100 mL bioreactor.

To visualize and perform the data treatment, download the Data_Treatment directory. After that open Matlab, add the path with subfolders in the Home tab and run the postProcessing script.

To perform your own simulations make sure you have OpenFOAM 8 installed. Then, download the reatorWALE cases and solver, and compile the solver by typping the wclean and wmake commands in the terminal (inside the solver directory). To set-up the mesh run the mesh/meshy bash file to create the mesh (in the case you want to change the opreational volume check point 4 first**). After select the operating conditions:

1 - Agitation rate (in rad/s) in constant/dynamicMeshDict;
2 - Kinematic viscosity and diffusivity as a function of temperature (which you can calculate, for water and oxygen, using the calcNuH2O and calcDab functions from the MATLAB directory) in constant/transportProperties;
3 - Saturation concentration as a function of temperature (which you can calculate, for water and oxygen, using the calcCsat function from the MATLAB directory) in 0/CO2 for the top patch;
**4 - Volume in the second column of the blockMeshDict vertices matrix (change only the largest value and redo the mesh. You can compare the values from the 60mL and 100mL and extrapolate a value based on the surface area to estimate the height you want).

You can also change the boundary conditions of the simulation, for example, the slip condition from the 0/U file by changing the slip coefficient (0 => slip; 1 => noSlip).
Additionally, you can change the refinement level of the mesh by changing the initial number of blocks in the blockMeshDict and changing several of the variables available in the hexMeshDict from the system folder (you can compare the reatorWALE and reatorWALE_refined blockMeshDict and hexMeshDict files for a general idea).

The mesh and meshy files can also be changed to include additional refinement near the surface.

To change the histogram parameters, start and end-time, as well as maxCFL access the system/controlDict. Here you can also add programmable functions such as the KMGS and EDR calculations.

To change the discretization schemes access system/fvSchemes.

To change the solution methods and parameters access system/fvSolution.

In case you run into any problems you can contact me through e-mail: pedro.miguel.neto@tecnico.ulisboa.pt
