function par = parameters_O2  
%Ref: Schmidt,R.;Wagner,W.: A new form of the equation of state for pure 
% substances and its application to oxygen. Fluid Phase Equilibirua.
% 1985, 19, 175. 
% The parameter tables were copied from: Leen van der Ham: An empirical 
%  Helmholtz energy based model for the calculation of thermodynamic 
%  properties of N2-Ar-O2 mixtures, NTNU, 2008

  par.species = 'O2';
  par.sections = [13,32];
  par.R    = 8314.46261815324;     % Universal gas constant (J/(kmol K))
  par.Mw   = 31.9988;       % Molar mass (kg/kmol)
  par.Tc   = 154.581;       % Critical temperature (K)
  par.rhoc = 13.63;        % Critical density (kmol/m3)
  par.vc   = 1/par.rhoc;    % Crtical molar volume (m3/kmol)
  par.pc   = 5.043e6;       % Critical pressure (Pa)
  par.Tt   = 54.361;        % Triple point temperature (K)
  par.pt   = 0.1463e5;      % Triple point pressure (Pa)


% van der ham Table 6
  par.ig_a = [
  -0.740775E-3
  -0.664930E-4
   0.250042E+1
  -0.214487E+2
   0.101258E+1
  -0.944365E+0
   0.145066E+2
   0.749148E+2
   0.414817E+1
];
  par.phi_ig = @(tau,delta,par,max_order)phi_ig_O2(tau,delta,par,max_order);
% van der ham Table 7
  A = [   1      3.983768749E-01    1       0
   2     -1.846157454E+00    1     1.5
   3      4.183473197E-01    1     2.5
   4      2.370620711E-02    2    -0.5
   5      9.771730573E-02    2     1.5
   6      3.017891294E-02    2       2
   7      2.273353212E-02    3       0
   8      1.357254086E-02    3       1
   9     -4.052698943E-02    3     2.5
  10      5.454628515E-04    6       0
  11      5.113182277E-04    7       2
  12      2.953466883E-07    7       5
  13     -8.687645072E-05    8       2
  14     -2.127082589E-01    1       5
  15      8.735941958E-02    1       6
  16      1.275509190E-01    2     3.5
  17     -9.067701064E-02    2     5.5
  18     -3.540084206E-02    3       3
  19     -3.623278059E-02    3       7
  20      1.327699290E-02    5       6
  21     -3.254111865E-04    6     8.5
  22     -8.313582932E-03    7       4
  23      2.124570559E-03    8     6.5
  24     -8.325206232E-04   10     5.5
  25     -2.626173276E-05    2      22
  26      2.599581482E-03    3      11
  27      9.984649663E-03    3      18
  28      2.199923153E-03    4      11
  29     -2.591350486E-02    4      23
  30     -1.259630848E-01    5      17
  31      1.478355637E-01    5      18
  32     -1.011251078E-02    5      23
  ];
  par.n = A(:,2);
  par.d = A(:,3);
  par.t = A(:,4);
  par.c = [ones(11,1)*2;ones(8,1)*4];
end

function res = phi_ig_O2(tau,delta,par,max_order)
% Ideal gas part of Helmholtz free energy (J/kmol)
% Source: Leen van der Ham, An empirical Helmholtz energy based model 
%   for the calculation of thermodynamic properties of N2-Ar-O2 mixtures
%   NTNU, December 20008
  a = par.ig_a;
  v0 = 24.44921;   % Reference molar volume at T0 = 298.15K, p0 = 101325pa
  delta0 = par.vc/v0;
  alpha = exp(a(7)*tau);
  beta7 = log(alpha-1);
  beta7_t = a(7)*alpha/(alpha-1);
  beta7_tt = -beta7_t^2/alpha;
 
  alpha8 = exp(a(8)*tau);
  beta8 = log(1+2/3/alpha8);;
  beta8_t = -2*a(8)/(3*alpha8+2);
  beta8_tt = 6*a(8)^2*alpha8/(3*alpha8+2)^2;

  phi = log(delta/delta0) + a(1)*tau^1.5 + a(2)*tau^-2 + a(3)*log(tau)...
    + a(4)*tau + a(5)*beta7 + a(6)*beta8;
  phi_t = 1.5*a(1)*tau^0.5 -2*a(2)*tau^-3 + a(3)/tau +a(4) +a(5)*beta7_t + a(6)*beta8_t;
  phi_d = 1/delta;
  phi_tt =  0.75*a(1)*tau^-0.5 +6*a(2)*tau^-4 -a(3)*tau^-2 ...
    +a(5)*beta7_tt + a(6)*beta8_tt;
  phi_dd = -delta^(-2);
  phi_td = 0;
  res = [phi,phi_t,phi_d,phi_tt,phi_td,phi_dd];
  if max_order > 2
    beta7_ttt = a(7)^3*alpha*(alpha+ 1)/(alpha-1)^3;
    beta8_ttt = -6*a(8)^3*alpha8*(3*alpha8-2)/(3*alpha8+2)^3;
    phi_ttt = -0.375*a(1)*tau^-1.5 -24*a(2)*tau^-5 +2*a(3)*tau^-3 ...
      + a(5)*beta7_ttt + a(6)*beta8_ttt;
    phi_ddd = 2*delta^(-3);
    res = [res,phi_ttt,0,0,phi_ddd];
  end
end

