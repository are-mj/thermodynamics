function par = CO2parameters
% Parameters for CO2 to be used with the CO2helmholtz.m
% and with the 'thermo' thermodynamic object.  
% Sources: 
%    General thermodynamic properties: NIST webbook
%    Thermodynaic model parameters:  
%      R. Span, W. Wagner,  J. Phys. Chem. Ref. Data 25, 1509 (1996);

  par.species = 'CO2';
  par.casno = '124389';
  par.sections = [7,34,39,42];
  par.R    = 8314.4626;     % Universl gas constant (J/(kmol K)
  par.Mw = 44.0098;         % Molar mass (kg/kmol)
  par.Tc   = 304.1282;      % Critical temperature (K)
  par.vc   = 0.094118;      % Crtical molar volume (NIST - m3/kmol)
  par.pc   = 73.77296e5;    % Critical pressure (Pa)
  par.rhoc = 1/par.vc;      % Critical density (kmol/m3)
  par.Tt = 216.592;         % Triple point temperature (K)
  par.pt = 5.1796e5;        % Triple point pressure
%   par.vlt = 0.037345;       % Liquid molar volume at triple point (m3/kmol)
%   par.vvt = 3.1982;         % Vapur phase molar volume at triple point
%   % Thermodynamic properties used by Span & Wagner:
%   par.R = 8314.51;
%   par.rhoc = 467.6/par.Mw;   
%   par.vc = 1/par.rhoc;
%   par.Tc = 304.1282;
%   par.pc = 7.3773e6;
%   par.Tt = 216.592;
%   par.pt = 5.1795e5;
  
  % Ideal gas (Table 27):
  par.ig_a = [
    8.37304456 
    -3.70454304
    2.50000000 
    1.99427042 
    0.62105248 
    0.41195293 
    1.04028922 
    0.08327678];
%   par.ig_a(1:2) = [-6.1248;15.9926];  % Change to NIST reference state
    par.ig_a(1:2) = par.ig_a(1:2) + [-120543;73334]/par.R;
  par.ig_b = [
    3.15163 
    6.11190 
    6.77708 
    11.32384
    27.08792];
  
  % Residual contribution  (Table 31)
  par.n = [
 0.38856823203161
 0.29385475942740e1
-0.55867188534934e1
-0.76753199592477
 0.31729005580416
 0.54803315897767 
 0.12279411220335
 0.21658961543220e1
 0.15841735109724e1
-0.23132705405503
 0.58116916431436e-1
-0.55369137205382 
 0.48946615909422
-0.24275739843501e-1 
 0.62494790501678e-1
-0.12175860225246
-0.37055685270086 
-0.16775879700426e-1 
-0.11960736637987 
-0.45619362508778e-1
 0.35612789270346e-1
-0.74427727132052e-2
-0.17395704902432e-2
-0.21810121289527e-1 
 0.24332166559236e-1 
-0.37440133423463e-1 
 0.14338715756878 
-0.13491969083286
-0.23151225053480e-1 
 0.12363125492901e-1 
 0.21058321972940e-2 
-0.33958519026368e-3 
 0.55993651771592e-2 
-0.30335118055646e-3 
-0.21365488688320e3 
 0.26641569149272e5 
-0.24027212204557e5 
-0.28341603423999e3 
 0.21247284400179e3 
-0.66642276540751
 0.72608632349897 
 0.55068668612842e-1];
 
par.d = [
1
1
1
1
2 
2 
3 
1
2 
4 
5 
5 
5
6 
6 
6 
1 
1 
4 
4 
4 
7 
8
2 
3 
3 
5
5
6
7 
8 
10 
4 
8 
2
2
2
3
3];

par.a = [
3.500 
3.500 
3.000];

par.t = [
0.00 
0.75 
1.00 
2.00 
0.75 
2.00 
0.75 
1.50 
1.50 
2.50 
0.00 
1.50 
2.00 
0.00 
1.00 
2.00 
3.00 
6.00 
3.00 
6.00 
8.00 
6.00 
0.00 
7.00 
12.00 
16.00 
22.00 
24.00 
16.00 
24.00 
8.00 
2.00 
28.00 
14.00
1.00 
0.00 
1.00 
3.00 
3.00 ];

par.b = [
0.875 
0.925 
0.875 ];

par.c = [ones(9,1);2*ones(7,1);3*ones(3,1);4*ones(6,1);5;6];

par.alpha = [
25 
25 
25 
15 
20];

par.beta = [
325 
300 
300 
275 
275 
0.300 
0.300 
0.300];

par.A = 0.7*ones(3,1);

par.gamma = [
1.16  
1.19  
1.19  
1.25  
1.22  ];
par.epsilon = ones(5,1);
par.B = [0.3;0.3;1];
par.C = [10;10;12.5];
par.D = 275*ones(3,1);

% Vapour-liquid saturation properties: 
% Pressure
par.as    = [-7.0602087, 1.9391218,-1.6463597,-3.2995634];
par.ase = [1;1.5;2;4]; 
% Liquid volume
par.bs = [1,1.83056,-1.55652,17.72784,-71.39499,170.18692,-237.26564,...
  178.57019,-55.65447];
par.bse = (0:8)'/3;
% Vapour volume
par.cs = [0.03520,  5.45521,-19.21051,140.49716,-564.97860,1473.50053,...
        -2300.80720,1981.73022,-695.9773];
par.cse = (0:8)'/2;

% Vapour/solid sublimation curve(eq. 3.12 in Span & Wagner)
par.sublima = [-14.740846, 2.4327015,  -5.3061778];
par.sublime = [1;1.9;2.9];  % Exponentials

% Liquid/solid melting curve
par.melta = [1955.5390 2055.4593];
par.melte = [1;2];  % Exponentials
end

  