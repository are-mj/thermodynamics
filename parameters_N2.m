function par = parameters_N2
% Ref: Span,R.;Lemmon,E.W.;Jacobsen,R.T.;Wagner,W.;Yokozeki,A.
%  A reference equation of state for the thermodynamic properties of  
%  nitrogen fortemperatures from 63.151 to 1000 K and pressures to 2200 MPa. 
%  J. Phys. Chem. Ref. Data.2000, 29, 1361. 
% The parameter tables were copied from: Leen van der Ham: An empirical 
%  Helmholtz energy based model for the calculation of thermodynamic 
%  properties of N2-Ar-O2 mixtures, NTNU, 2008

  par.species = 'N2';
  par.sections = [6,32,36];
  par.R    = 8314.46261815324;     % Universal gas constant (J/(kmol K))
  par.Mw   = 28.0135;       % Molar mass (kg/kmol)
  par.Tc   = 126.192;       % Critical temperature (K)
  par.rhoc = 11.184;        % Critical density (kmol/m3)
  par.vc   = 1/par.rhoc;    % Crtical molar volume (m3/kmol)
  par.pc   = 3.396e6;       % Critical pressure (Pa)
  par.Tt   = 63.151;        % Triple point temperature (K)
  par.pt   = 12.523e5;      % Triple point pressure (Pa)


% van der ham Table 6
  par.ig_a = [
  2.5E+0
  12.76952708E+0
   -0.784163E-2
    -1.934819E-4
    -1.247742E-5
     6.678326E-8
     1.012941E+0
     26.65788E+0];
  par.phi_ig = @(tau,delta,par,max_order)phi_ig_N2(tau,delta,par,max_order);
  
% van der ham Table 7
  A = [   1    9.24803575275E-01    1    0.25   0
     2   -4.92448489428E-01    1   0.875   0
     3    6.61883336938E-01    2     0.5   0
     4   -1.92902649201E+00    2   0.875   0
     5   -6.22469309629E-02    3   0.375   0
     6    3.49943957581E-01    3    0.75   0
     7    5.64857472498E-01    1     0.5   1
     8   -1.61720005987E+00    1    0.75   1
     9   -4.81395031883E-01    1       2   1
    10    4.21150636384E-01    3    1.25   1
    11   -1.61962230825E-02    3     3.5   1
    12    1.72100994165E-01    4       1   1
    13    7.35448924933E-03    6     0.5   1
    14    1.68077305479E-02    6       3   1
    15   -1.07626664179E-03    7       0   1
    16   -1.37318088513E-02    7    2.75   1
    17    6.35466899859E-04    8    0.75   1
    18    3.04432279419E-03    8     2.5   1
    19   -4.35762336045E-02    1       4   2
    20   -7.23174889316E-02    2       6   2
    21    3.89644315272E-02    3       6   2
    22   -2.12201363910E-02    4       3   2
    23    4.08822981509E-03    5       3   2
    24   -5.51990017984E-05    8       6   2
    25   -4.62016716479E-02    4      16   3
    26   -3.00311716011E-03    5      11   3
    27    3.68825891208E-02    5      15   3
    28   -2.55856846220E-03    8      12   3
    29    8.96915264558E-03    3      12   4
    30   -4.41513370350E-03    5       7   4
    31    1.33722924858E-03    6       4   4
    32    2.64832491957E-04    9      16   4
    33    1.96688194015E+01    1       0   2
    34   -2.09115600730E+01    1       1   2
    35    1.67788306989E-02    3       2   2
    36    2.62767566274E+03    2       3   2];
  par.n = A(:,2);
  par.d = A(:,3);
  par.t = A(:,4);
  par.c = A(par.sections(1)+1:par.sections(2),5);

  par.alpha = [20,20,15,25]';
  par.beta = [325,325,300,275]';
  par.gamma = [1.16,1.16,1.13,1.25]';
  par.epsilon = [1,1,1,1]';
  
end

function res = phi_ig_N2(tau,delta,par,max_order)
% Ideal gas part of Helmholtz free energy (J/kmol)
% Source: Leen van der Ham, An empirical Helmholtz energy based model 
%   for the calculation of thermodynamic properties of N2-Ar-O2 mixtures
%   NTNU, December 20008
  a = par.ig_a;
  alpha = exp(a(8)*tau);
  beta = log(1-1/alpha);
  beta_t = a(8)/(alpha-1);
  beta_tt = -beta_t^2*alpha;
  
  phi = log(delta) + a(1)*log(tau) + a(2) + a(3)*tau  + a(4)*tau^-1 ...
    + a(5)*tau^-2 + a(6)*tau^-3 + a(7)*beta;
  phi_t = a(1)*tau^-1 + a(3) -a(4)*tau^-2 -2*a(5)*tau^-3 -3*a(6)*tau^-4 ...
    +a(7)*beta_t; 
  phi_d = 1/delta;
  phi_tt =  -a(1)*tau^-2 +2*a(4)*tau^-3 +6*a(5)*tau^-4 +12*a(6)*tau^-5 ...
    +a(7)*beta_tt;
  phi_dd = -delta^(-2);
  phi_td = 0;
  res = [phi,phi_t,phi_d,phi_tt,phi_td,phi_dd];
  if max_order > 2
    beta_ttt = beta_t^3*alpha*(alpha+1);
    phi_ttt = 2*a(1)*tau^-3 -6*a(4)*tau^-4 -24*a(5)*tau^-5 ...
      -60*a(6)*tau^-6 + a(7)*beta_ttt;
    phi_ddd = 2*delta^(-3);
    res = [res,phi_ttt,0,0,phi_ddd];
  end
end

