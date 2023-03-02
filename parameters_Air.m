function par = parameters_Air  
%Parameters for termodynamic model for air as a pseudo species
%Ref: Lemmon,  E.W.;  Jacobsen,  R.T.;  Penoncello,  S.G.;  Friend,  D.G. 
% Thermodynamic  properties  of  air  and  mixtures  of  nitrogen,  argon, 
% and oxygen from 60 to 2000 K at pressures to 2000 MPa. 
% J. Phys. Chem. Ref. Data. 2000, 29, 331. 
%The parameter tables were copied from: Leen van der Ham: An empirical 
%  Helmholtz energy based model for the calculation of thermodynamic 
%  properties of N2-Ar-O2 mixtures, NTNU, 2008

  par.species = 'Air';
  x_air = [0.7812,0.2096,0.0092];  % N2, O2 and Ar mole fractions
  Mw = [28.0135,31.9988,39.9480];  % N2, O2 and Ar molar masses
  par.sections = [10,19];
  par.R    = 8314.46261815324;     % Universal gas constant (J/(kmol K))
  par.Mw   = Mw*x_air';            % Molar mass (kg/kmol)
  % Note that the maxcondentherm replaces the critical point in teh model
  % Maxcondentherm: the maximum temperature point of the phase envenlope 
  par.Tc   = 132.6312;      % Maxcondentherm (K)
  par.rhoc = 10.4477;       % Maxcondentherm density (kmol/m3)
  par.vc   = 1/par.rhoc;    % Maxcondentherm (m3/kmol)
  par.pc   = 3.78502e6;     % Maxcondentherm pressure (Pa)
  par.Tt   = 63.151;        % Triple point temperature (K)
  par.pt   = 12.523e5;      % Triple point pressure (Pa)


% Lemmon table 12
  par.ig_a = [ 0.605719400e-7
-0.210274769e-4
-0.158860716e-3
-13.841928076
 17.275266575
-0.195363420e-3
 2.490888032
 0.791309509
 0.212236768
-0.197938904
 25.36365
 16.90741
 87.31279
];

par.phi_ig = @(tau,delta,par,max_order)phi_ig_Air(tau,delta,par,max_order);
% Lemmon table 13:
A = [
   1       0.118160747229    1       0   0
   2       0.713116392079    1    0.33   0
   3   -0.161824192067e+1    1    1.01   0
   4    0.714140178971e-1    2       0   0
   5   -0.865421396646e-1    3       0   0
   6       0.134211176704    3    0.15   0
   7    0.112626704218e-1    4       0   0
   8   -0.420533228842e-1    4     0.2   0
   9    0.349008431982e-1    4    0.35   0
  10    0.164957183186e-3    6    1.35   0
  11      -0.101365037912    1     1.6   1
  12      -0.173813690970    3     0.8   1
  13   -0.472103183731e-1    5    0.95   1
  14   -0.122523554253e-1    6    1.25   1
  15      -0.146629609713    1     3.6   2
  16   -0.316055879821e-1    3       6   2
  17    0.233594806142e-3   11    3.25   2
  18    0.148287891978e-1    1     3.5   3
  19   -0.938782884667e-2    3      15   3
  ];

  par.n = A(:,2);
  par.d = A(:,3);
  par.t = A(:,4);
  par.c = A(par.sections(1)+1:par.sections(2),5);
end

function res = phi_ig_Air(tau,delta,par,max_order)
% Ideal gas part of Helmholtz free energy (J/kmol)
% Ref: Lemmon,  E.W.;  Jacobsen,  R.T.;  Penoncello,  S.G.;  Friend,  D.G. 
% Thermodynamic  properties  of  air  and  mixtures  of  nitrogen,  argon, 
% and oxygen from 60 to 2000 K at pressures to 2000 MPa. 
% J. Phys. Chem. Ref. Data. 2000, 29, 331. 
% Equation 24 and Table 12
  N = par.ig_a;
  alpha = exp(-N(11:13)*tau);
  phi = log(delta) + tau.^[-3:1,1.5]*N(1:6)+N(7)*log(tau) ...
    + N(8:9)'*log(1-alpha(1:2)) + N(10)*log(2/3+alpha(3)^-1);
  phi_t = [-3:1,1.5].*N(1:6)'*tau.^([-3:1,1.5]-1)' + N(7)*tau^-1 ...
    + N(8)*N(11)*alpha(1)/(1-alpha(1)) ...
    + N(9)*N(12)*alpha(2)/(1-alpha(2)) ...
    + N(10)*N(13)/(2/3*exp(-N(13)*tau)+1);
  phi_tt = [12,6,2].*N(1:3)'*tau.^[-5,-4,-3]' ...
    + 0.75*N(6)*tau^-0.5 -N(7)*tau^-2 ...
    - N(8)*N(11)^2*alpha(1)/(1-alpha(1))^2 ...
    - N(9)*N(12)^2*alpha(2)/(1-alpha(2))^2 ...
    + N(10)*N(13)^2*6*alpha(3)/(3+2*alpha(3))^2;
  phi_d = 1/delta;
  phi_dd = -1/delta^2;
  res = [phi,phi_t,phi_d,phi_tt,0,phi_dd];
  if max_order > 2
    phi_ttt = -[60 24 6].*N(1:3)'*tau.^[-6,-5,-4]' ...
    - 0.375*N(6)*tau^-1.5 +2*N(7)*tau^-3 ...
    + N(8)*N(11)^3*alpha(1)*(1+alpha(1))/(1-alpha(1))^3 ...
    + N(9)*N(12)^3*alpha(2)*(1+alpha(2))/(1-alpha(2))^3 ...
    - N(10)*N(13)^3*6*alpha(3)*(3-2*alpha(3))/(3+2*alpha(3))^3;
    phi_ddd = 2*delta^-3;
    res = [res,phi_ttt,0,0,phi_ddd];
  end
end
