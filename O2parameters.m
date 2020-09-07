function par = N2parameters  
% Parameters from van der Ham
  par.species = 'O2';
  par.sections = [13,32];
  par.R    = 8314.3714;     % Universl gas constant (J/(kmol K)
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




