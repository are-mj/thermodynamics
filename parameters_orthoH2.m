function par = parameters_orthoH2()
% Parameters for ortho H2

% General parameters, taken from Leachman et al.:
  par.species = 'orthoH2';
  par.casno = '';
  par.R    = 8314.462618;   % Universl gas constant (J/(kmol K)
  par.Tc   = 33.22;         % Critical temperature (K)
  par.pc   = 13.1065e5;     % Critical pressue (Pa)
  par.rhoc = 15.445;        % Crititcal molar density (kmol/m3) 
  par.vc   = 1/par.rhoc;    % Crititcal molar volume (m3/kmol)
  par.Mw   = 2.01588;       % Molar mass (kg/kmol)
  par.Tt   = 14.008;        % Triple point temperature, K
  par.pt   = 7461;          % Triple point pessure (K)

% Parameters for the Helmholtz equation of state for normal H2 from 
% J. W. Leachman & al., "Fundamental Equations of State for Parahydrogen,
%  Normal Hydrogen, and Orthohydrogen", Journal of Physical and Chemical 
%  Reference Data · September 2009, pp. 721-784

  % Table 4:
  par.ig_a = [-1.4675442336 1.8845068862 1.5 2.54151 -2.3661 1.00365 ...
    1.22447]';
  par.ig_b = [25.7676098736 43.4677904877 66.0445514750 209.7531607465]';
  
  par.sections = [7,9,14];
  % Table 5:
  par.n = [-6.83148 0.01 2.11505 4.38353 0.211292 -1.00939 0.142086 ...
    -0.87696 0.804927 -0.710775 0.0639688 0.0710858 -0.087654 0.647088]';
  par.t = [0.7333 1 1.1372 0.5136 0.5638 1.6248 1.829 2.404 2.105 4.1 ...
    7.658 1.259 7.589 3.946]';
  par.d = [1 4 1 1 2 2 3 1 3 2 1 3 1 1]';
  par.c = [1 1]';
  
  % Table 6:
  par.alpha = [1.169;0.894;0.04;2.072;1.306];
  par.beta = [0.4555;0.4046;0.0869;0.4415;0.5743];
  par.gamma = [1.5444;0.6627;0.763;0.6587;1.4327];
  par.epsilon = [0.6366;0.3876;0.9437;0.3976;0.9626];
  
  % Saturation pressure data (Table 8)   (ps = pc*exp(Tc/T*as*theta.^ase);  
  par.as = [-4.88684 1.05310 0.856947 -0.185355];
  par.ase = [1 1.5 2.7 6.2]';
  % Saturated liquid molar volume;
  par.bs = [1.001220,3.295132,-4.129730,4.980071,-2.504104];
  par.bse = (0:4)'/2;
  % Saturated vapour molar volume;
  par.cs = [-2.898846,51.893600,-376.027881,1527.104902,-3587.700043,4914.542445,-3638.389436,1132.737750];
  par.cse = (0:7)'/3;
  
  par.cs = [0.000591,2.986426,6.038556,-43.811536,149.135328,-210.019433,116.302058];
  par.cse = (0:6)'/2;
end
