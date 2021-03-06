function par = parameters_Ar
%Ref: Tegeler,Ch.;Span,R.;Wagner,W.: A new equation of state for argon 
% covering the fluid region for temperatures from the melting line to 
% 700 K at pressures up to 1000 MPa. 
% J. Phys. Chem. Ref. Data. 1999, 28, 779. 
% The parameter tables were copied from: Leen van der Ham: An empirical 
%  Helmholtz energy based model for the calculation of thermodynamic 
%  properties of N2-Ar-O2 mixtures, NTNU, 2008

  par.species = 'Ar';
  par.sections = [12,37,41];
  par.R    = 8314.46261815324;     % Universal gas constant (J/(kmol K))
  par.Mw   = 39.9480;       % Molar mass (kg/kmol)
  par.Tc   = 150.687;       % Critical temperature (K)
  par.rhoc = 13.407;        % Critical density (kmol/m3)
  par.vc   = 1/par.rhoc;    % Crtical molar volume (m3/kmol)
  par.pc   = 4.863e6;       % Critical pressure (Pa)
  par.Tt   = 83.806;        % Triple point temperature (K)
  par.pt   = 68.891e5;      % Triple point pressure (Pa)


% van der ham Table 6
  par.ig_a = [
   8.31666243
  -4.94651164
   1.5
];

% van der ham Table 7
  A = [   1  8.8722304990011E-02    1       0   0
   2  7.0514805167298E-01    1    0.25   0
   3 -1.6820115654090E+00    1       1   0
   4 -1.4909014431486E-01    1    2.75   0
   5 -1.2024804600940E-01    1       4   0
   6 -1.2164978798599E-01    2       0   0
   7  4.0035933626752E-01    2    0.25   0
   8 -2.7136062699129E-01    2    0.75   0
   9  2.4211924579645E-01    2    2.75   0
  10  5.7889583185570E-03    3       0   0
  11 -4.1097335615341E-02    3       2   0
  12  2.4710761541614E-02    4    0.75   0
  13 -3.2181391750702E-01    1       3   1
  14  3.3230017695794E-01    1     3.5   1
  15  3.1019986287345E-02    3       1   1
  16 -3.0777086002437E-02    4       2   1
  17  9.3891137419581E-02    4       4   1
  18 -9.0643210682031E-02    5       3   1
  19 -4.5778349276654E-04    7       0   1
  20 -8.2659729025197E-05   10     0.5   1
  21  1.3013415603147E-04   10       1   1
  22 -1.1397840001996E-02    2       1   2
  23 -2.4455169960535E-02    2       7   2
  24 -6.4324067175955E-02    4       5   2
  25  5.8889471093674E-02    4       6   2
  26 -6.4933552112965E-04    8       6   2
  27 -1.3889862158435E-02    3      10   3
  28  4.0489839296910E-01    5      13   3
  29 -3.8612519594749E-01    5      14   3
  30 -1.8817142332233E-01    6      11   3
  31  1.5977647596482E-01    6      14   3
  32  5.3985518513856E-02    7       8   3
  33 -2.8953417958014E-02    7      14   3
  34 -1.3025413381384E-02    8       6   3
  35  2.8948696775778E-03    9       7   3
  36 -2.2647134304796E-03    5      24   3
  37  1.7616456196368E-03    6      22   4
  38  5.8552454482774E-03    2       3   0
  39 -6.9251908270028E-01    1       1   0
  40  1.5315490030516E+00    2       0   0
  41 -2.7380447449783E-03    3       0   0
  ];
  par.n = A(:,2);
  par.d = A(:,3);
  par.t = A(:,4);
  par.c = A(13:37,5);
  
  par.alpha = [20 20 20 20]';
  par.beta = [250 375 300 225]';
  par.gamma = [1.11 1.14 1.17 1.11]';
  par.epsilon = ones(4,1);
end




