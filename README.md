# thermodynamics

This repository contains Matlab implementations of high-accuracy thermodynamic models from the litterature, together with the 'thermo' class for calculating thermodynamic variables.  Currenty, available species are H2, paraH2, orthoH2, N2, O2, Ar, H2O, and CO2, together with Air treated as a single pseudo species. 

**Main files:**

thermo.m:  Thermodynamic object that contains methods and properties that enable the calculation of thermodynamic variables and processes.

helmholtz.m:  Helmholtz molar free energy and partial derivatives.  Used by thermo.Tvcalc

parametersXX.m: Parameters for spescies XX.  Used by thermo and helmholtz.

A number of supplementary funcrions and examples are also included.

With the exception of Air and orthoH2, the models and parameters are the same as used by NIST for calculating Thermophysical Properties of Fluid Systems
   https://webbook.nist.gov/chemistry/fluid/

Use the Matlab help command for details (e.g. help thermo)

**Documentation**:

- User guide for Matlab class thermo.pdf

- Properties from Helmholtz.pdf

- Shock tube model for real gases.pdf


[![View H2 and CO2 thermodynamics on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/73950)

