function par = parameters_H2
% Parameters for normal H2: 75% ortho-H2, 25% para-G2
% This is the high temperature equlibrium composition
% (ortho-H2 > 74.8% for T > 267K)

% General parameters, taken from NIST:
  par.species = 'H2';
  par.casno = '1333740';
  par.R    = 8314.46261815324;     % Universal gas constant (J/(kmol K))
  par.Tc   = 33.145;        % Critical temperature (K)
  par.pc   = 12.964e5;      % Critical pressue (Pa)
  par.rhoc = 15.508;        % Crititcal molar density (kmol/m3) 
  par.vc   = 1/par.rhoc;    % Crititcal molar volume (m3/kmol)
  par.Mw   = 2.01588;       % Molar mass (kg/kmol)
  par.Tt   = 13.95;         % Triple point temperature, K
  par.pt   = 7210;          % Triple point pessure (K)

% Parameters for the Helmholtz equation of state for normal H2 from 
% J. W. Leachman & al., "Fundamental Equations of State for Parahydrogen,
%  Normal Hydrogen, and Orthohydrogen", Journal of Physical and Chemical 
%  Reference Data · September 2009, pp. 721-784

  % Table 4:
  par.ig_a = [-1.4579856475,1.888076782,1.5,1.616 ,-0.4117 ,-0.792 ,0.758 , ...
    1.217]';
  par.ig_b = [16.0205159149,22.6580178006,60.0090511389,...
    74.9434303817,206.9392065168]';

  par.sections = [7,9,14];
  % Data from Leachman's tables 5 and 6 as copied from Table 2 in Yang: 
  % A thermodynamic analysis of refueling of a hydrogen tank
  % Int. J. Hydrogen Energy, 34 (2009) 6712–6721
  nn = 0;  % Missing values
  A =[
  %   N           t       d c   phi  beta  gamma      D
  1  -6.93643     0.6844  1 0   nn     nn     nn     nn
  2   0.01        1       4 0   nn     nn     nn     nn
  3   2.1101      0.989   1 0   nn     nn     nn     nn
  4   4.52059     0.489   1 0   nn     nn     nn     nn
  5   0.732564    0.803   2 0   nn     nn     nn     nn
  6  -1.34086     1.1444  2 0   nn     nn     nn     nn
  7   0.130985    1.409   3 0   nn     nn     nn     nn
  8  -0.777414    1.754   1 1   nn     nn     nn     nn
  9   0.351944    1.311   3 1   nn     nn     nn     nn
  10 -0.0211716   4.187   2 nn -1.685 -0.171  0.7164 1.506
  11  0.0226312   5.646   1 nn -0.489 -0.2245 1.3444 0.156
  12  0.032187    0.791   3 nn -0.103 -0.1304 1.4517 1.736
  13 -0.0231752   7.249   1 nn -2.506 -0.2785 0.7204 0.67
  14  0.0557346   2.986   1 nn -1.607 -0.3967 1.5445 1.662];
  par.n = A(:,2);
  par.t = A(:,3);
  par.d = A(:,4);
  par.c = A(8:9,5);
  par.alpha = -A(10:14,6);
  par.beta =  -A(10:14,7);
  par.gamma =  A(10:14,8);
  par.epsilon =  A(10:14,9);

  % Saturation pressure data (Table 8)   (ps = pc*exp(Tc/T*as*theta.^ase);  
  par.as = [-4.89789, 0.988558, 0.349689, 0.499356];
  par.ase = [1;1.5;2;2.85];
  % Saturated liquid volume  (vl = vc/(bs*theta.^bse)
  par.bs = [1,0.27768,12.94487,-61.00770,170.61820,-290.47816,298.08997,...
        -169.80272, 41.11494];
  par.bse = (0:8)'/3;
  % Saturated liquid volume (vv = vc/(cs*theta.^cse)
  par.cs = [0.014155,3.14575,9.14921,-128.63980,883.85951,-3521.00881,...
    9002.68518,-14926.18150,15586.28934,-9337.23784,2462.54301];
  par.cse = (0:10)'/2;