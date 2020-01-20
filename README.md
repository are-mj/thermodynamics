# thermodynamics

This repository contains Matlab implementations of high-accuracy thermodynamic models from the litterature.  So far for hydrogen and CO2. I also supply a Matlab thermodynamic object 'thermo' that simplifies working with the models.

Main files:

CO2helmholtz.m: Helmholtz molar free energy and partial derivatives for pure CO2.  
   Ref.: Span and Wagner, J. Phys. Chem. Reference Data 25, 1509 (1996)
   https://doi.org/10.1063/1.555991
   
H2helmholtz.m: Helmholtz molar free energy and partial derivatives for pure hydrogen.
   Ref: Leachman, J.W. et al., J. Phys. Chem. Reference Data 38, 721 (2009)

Both models are the same as used by NIST for calculating Thermophysical Properties of Fluid Systems
   https://webbook.nist.gov/chemistry/fluid/

thermo.m:  Thermodynamic object that contains methods and properties that enable the calculation of thermodynamic variables and processes.

Use the Matlab help command for details (e.g. help CO2helmholtz)

Thanks to PS for the InterX function used in H2twophase_region and CO2twophase_region. https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections
