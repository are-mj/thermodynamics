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