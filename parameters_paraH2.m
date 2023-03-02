function par = parameters_paraH2
% Parameters for para H2

% General parameters, taken from NIST and Leachman et al.:
  par.species = 'paraH2';
  par.casno = '';
  par.R    = 8314.462618;   % Universl gas constant (J/(kmol K)
  par.Tc   = 32.938;        % Critical temperature (K)
  par.pc   = 12.858e5;      % Critical pressue (Pa)
%   par.pc   = 12.8377e5;     % Critical pressue (Pa)  (NIST)
  par.rhoc = 15.538;        % Crititcal molar density (kmol/m3) 
  par.vc   = 1/par.rhoc;    % Crititcal molar volume (m3/kmol)
  par.Mw   = 2.01588;       % Molar mass (kg/kmol)
  par.Tt   = 13.8033;       % Triple point temperature, K
  par.pt   = 7041;          % Triple point pessure (K)

% Parameters for the Helmholtz equation of state for normal H2 from 
% J. W. Leachman & al., "Fundamental Equations of State for Parahydrogen,
%  Normal Hydrogen, and Orthohydrogen", Journal of Physical and Chemical 
%  Reference Data · September 2009, pp. 721-784

  % Table 4:
  par.ig_a = [-1.4485891134 1.884521239 1.5 4.30256 13.0289 -47.7365 ...
    50.0013 -18.6261 0.993973 0.536078]';
  par.ig_b = [15.1496751472 25.0925982148 29.4735563787 35.4059141417 ... 
    40.724998482 163.7925799988 309.2173173842 ]';
  
  par.sections = [7,9,14];
  % Table 5:
  par.n = [-7.33375 0.01 2.60375 4.66279 0.682390 -1.47078 0.135801 ...
    -1.05327 0.328239 -0.0577833 0.0449743 0.0703464 -0.0401766 0.119510]';
  par.t = [0.6855 1 1 0.489 0.774 1.133 1.386 1.619 1.162 3.96 5.276 ...
    0.99 6.791 3.19]';
  par.d = [1 4 1 1 2 2 3 1 3 2 1 3 1 1]';
  par.c = [1 1]';
  
  % Table 6
  par.alpha = [1.7437;0.5516;0.0634;2.1341;1.777];
  par.beta = [0.194;0.2019;0.0301;0.2383;0.3253];
  par.gamma = [0.8048;1.5248;0.6648;0.6832;1.493];
  par.epsilon = [1.5487;0.1785;1.28;0.6319;1.7104];
  
  % Saturation pressure data (Table 8)   (ps = pc*exp(Tc/T*as*theta.^ase);  
  par.as = [-4.87767 1.03359 0.82668 -0.129412];
  par.ase = [1 1.5 2.65 7.4]';
  
  % Saturated liquid molar volume;
  par.bs = [1.000253,1.215997,1.201567,-0.675420];
  par.bse = (0:3)'/3;
  % Saturated vapour molar volume;
%   par.cs = [1.000008,-0.918816,-4.273750,13.796604,-23.935002,22.334254,-8.000405];
  par.cse = (0:6)'/2;
  par.cs = [0.000000,2.685390,12.266407,-76.276128,223.099876,-288.454916,148.036040,];
end