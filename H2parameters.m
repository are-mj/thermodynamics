function par = H2parameters()
% Parameters for normal H2

% General parameters, taken from NIST:
  par.species = 'H2';
  par.R    = 8314.462618;   % Universl gas constant (J/(kmol K)
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
  par.a = [-1.4579856475 ,1.888076782 ,1.616 ,-0.4117 ,-0.792 ,0.758 , ...
    1.217]';
  par.b = [0,0,-16.0205159149,-22.6580178006,-60.0090511389,...
    -74.9434303817,-206.9392065168]';

  % Data from Leachman's tables 5 and 6 as copied from Table 2 in Yang: 
  % A thermodynamic analysis of refueling of a hydrogen tank
  % Int. J. Hydrogen Energy, 34 (2009) 6712–6721
  nn = 0;  % Missing values
  A =[
  %   N           t       d p   phi  beta  gamma      D
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
  par.N = A(:,2);
  par.t = A(:,3);
  par.d = A(:,4);
  par.phi = A(10:14,6);
  par.beta =  A(10:14,7);
  par.gamma =  A(10:14,8);
  par.D =  A(10:14,9);

  % Saturation pressure data (Table 8):
  par.Ns = [-4.89789; 0.988558; 0.349689; 0.499356];
  par.ks = [1;1.5;2;2.85];
  % Saturation pressure polynomial:
  par.pp = [
     1.784921538507186
    -5.778743203156100
     7.983636496366720
    -6.239109613294445
     3.035909635915719
    -0.968009110994933
     0.200502034945211
    -0.032242453027793
    -0.000225217801689
    -0.004831078987840
                     0]*1e3;
                   
  % Approximate saturated liquid volume:  
  %  vl = [theta-0.01,1,theta,thets^2,theta^3]*pvl   
  %  theta = 1 - T/Tc
   par.pvl = [
   0.000201307620658
   0.040037638985877
  -0.060940347137854
   0.105028421110417
  -0.073516224233979];     

% Approximate saturated vapour volume
%  vv = exp([theta-0.01,1,theta,thets^2,theta^3]*pvv)
  par.pvv = [
    -0.003852330905278
    -2.276541662489586
     7.729823951843128
    -8.822575591555243
    17.926208944187930];
end