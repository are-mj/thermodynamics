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