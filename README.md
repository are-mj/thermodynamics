# thermodynamics

This repository contains Matlab implementations of high-accuracy thermodynamic models from the litterature, together with the 'thermo' class for calculating thermodynamic variables.  Currenty, available species are H2, N2, O2, Ar, H2O, and CO2, together with Air treated as a single pseudo species. 

Main files:

thermo.m:  Thermodynamic object that contains methods and properties that enable the calculation of thermodynamic variables and processes.

helmholtz.m:  Helmholtz molar free energy and partial derivatives.  Used by thremo.Tvcalc

CO2parameters.m: Thermodynamic and model parameters for pure CO2.  
   Ref.: Span and Wagner, J. Phys. Chem. Reference Data 25, 1509 (1996)
   https://doi.org/10.1063/1.555991
   
H2parameters.m: Thermodynamic and model parameters for pure hydrogen.
   Ref: Leachman, J.W. et al., J. Phys. Chem. Reference Data 38, 721 (2009)
   
H2Oparameters.m: Thermodynamic and model parameters for pure H2O.
   Ref: W. Wagner, A. Pruss. Journal of Physical and Chemical Reference Data, 2002: 387-535. 

N2parameters.m, O2parameters.m, Arparameters.m: Thermodynamic and model parameters for pure nitrogen, oxygen and argon
   Ref: Lemmon et al.:  Phys. Chem. Reference Data 29, 331 (2000)
   
Airparameters.m: Thermodynamic and model parameters for air as a singke pseudo species
   Ref: Lemmon et al.:  Phys. Chem. Reference Data 29, 331 (2000)

With the exception of Air, the model and parameters are the same as used by NIST for calculating Thermophysical Properties of Fluid Systems
   https://webbook.nist.gov/chemistry/fluid/

Use the Matlab help command for details (e.g. help thermo)


[![View H2 and CO2 thermodynamics on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/73950)

