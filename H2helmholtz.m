function res = H2helmholtz(T,v,max_order)
% Helmholtz molar free energy and partial derivatives for pure hydrogen
%  T:   Temperature (K)
%  v:   Molar volume (m3/kmol)
%  par: Model paramters struct
%  max_order: Highest order of partial derivatives. Default: 2, Max: 3
% Model from Leachman, J.W. et al., 
%       J. Physical and Chemical Reference Data,vol. 38, p721 (2009)

  if nargin < 3
    max_order = 2;
  end
  par = H2parameters;

  % Reduced temperature and density:
  tau = par.Tc/T;
  delta = par.vc/v;
  tau_T = -par.Tc/T^2;
  tau_TT = 2*par.Tc/T^3;
  delta_v = -1/(v^2*par.rhoc);
  delta_vv = 2/(v^3*par.rhoc);
  delta_vvv = -6/(v^4*par.rhoc);

  % Ideal gas contribution
  res_ig = a_ig(tau,delta,par,max_order);

  % Residual contribution
  res_r = ar(tau,delta,par,max_order);

  res = res_ig + res_r;

  % reducet helmholts free energy and deriventive w.r.t tau (t) and delta (d)
  a = res(1);
  a_t = res(2);
  a_d = res(3);
  a_tt = res(4);
  a_td = res(5);
  a_dd = res(6);
  if max_order > 2
    a_ttt = res(7);
    a_ttd = res(8);
    a_tdd = res(9);
    a_ddd = res(10);
  end

  R = par.R;
  f    = R*T*a;
  f_T = R*(a + T*a_t*tau_T);
  f_v = R*T*a_d*delta_v;
  % f_TT = R*(2*a_t*tau_T + T*a_tt*tau_T^2+T*a_t*tau_TT);
  %  But terms 1 and 3 cancel, so:
  f_TT = R*T*a_tt*tau_T^2;
  f_Tv = R*(a_d*delta_v + T*a_td*tau_T*delta_v);
  f_vv = R*T*(a_dd*delta_v^2 + a_d*delta_vv);
  res = [f,f_T,f_v,f_TT,f_Tv,f_vv];
  if max_order > 2
    f_TTT = R*(a_tt*tau_T^2 + T*a_ttt*tau_T^3 + 2*T*a_tt*tau_T*tau_TT);
    f_TTv = R*T*a_ttd*tau_T^2*delta_v;
    f_Tvv = R*(a_dd*delta_v^2 + T*a_tdd*delta_v^2*tau_T + a_d*delta_vv ...
      + T*a_td*delta_vv*tau_T);
    f_vvv = R*(T*a_ddd*delta_v^3 + 3*T*a_dd*delta_v*delta_vv ...
      + T*a_d*delta_vvv);
    res = [res,f_TTT,f_TTv,f_Tvv,f_vvv];
  end
end

function res = a_ig(tau,delta,par,max_order)
% Reduced ideal gas Helmholtz free energy with derivatives
% Ref.: Leachman, J.W. et al., 
%       J. Physical and Chemical Reference Data,vol. 38, p721 (2009)
  A = par.a(3:7);
  B = par.b(3:7);
  y = exp(B*tau);
	a = log(delta)+1.5*log(tau)+par.a(1)+par.a(2)*tau + A'*log(1-y);
  a_T = 1.5/tau +par.a(2) - (A.*B)'*(y./(1-y));
  a_d = 1/delta;
  a_TT = -1.5/tau^2 - (A.*B.^2)'*(y./(1-y).^2);
  a_dd = -delta^(-2);
  a_Td = 0;
  res = [a,a_T,a_d,a_TT,a_Td,a_dd];
  if max_order > 2
    a_TTT = 3/tau^3-(A.*B.^3)'*(y.*(1+y)./(1-y).^3);
    a_ddd = 2*delta^(-3);
    res = [res,a_TTT,0,0,a_ddd];
  end
end

function res = ar(tau,delta,par,max_order)
% Reduced residual Helmholtz free energy a with derivatives
% Ref.: Leachman, J.W. et al., 
%       J. Physical and Chemical Reference Data,vol. 38, p721 (2009)
% res = [a,a_t,a_d,a_tt,a_td,a_dd,a_ttt,a_ttd,a_tdd,a_ddd]

  N = par.N;
  t = par.t;
  d = par.d;
  phi = par.phi;
  beta = par.beta;
  gamma = par.gamma;
  D = par.D;
  
  s = N.*delta.^d.*tau.^t;
  s_t = s.*t/tau;
  s_d = s.*d/delta;
  s_tt = s_t.*(t-1)/tau;
  s_td = s_t.*d/delta;
  s_dd = s_d.*(d-1)/delta; 
  
  i1 = 1:7;
  S1 = sum(s(i1));
  S1_t = sum(s_t(i1));
  S1_tt = sum(s_tt(i1));
  S1_d = sum(s_d(i1));
  S1_dd = sum(s_dd(i1));
  S1_td = sum(s_td(i1));
  
  res1 = [S1,S1_t,S1_d,S1_tt,S1_td,S1_dd];
  
  % e = exp(-delta.^p(i2)), but p(i2) = 1:
  e = exp(-delta);
  i2 = 8:9;
  S2 = sum(s(i2))*e;
  S2_t = sum(s_t(i2))*e;
  S2_d = sum(s_d(i2)-s(i2))*e;
  S2_tt = sum(s_tt(i2))*e;
  S2_td = sum(t(i2)/tau.*(s_d(i2)-s(i2)))*e;
  S2_dd = sum(s_dd(i2)-2*s_d(i2)+s(i2))*e;
  
  res2 = [S2,S2_t,S2_d,S2_tt,S2_td,S2_dd];
 
  o = exp(phi.*(delta-D).^2+beta.*(tau-gamma).^2);
  o_t = 2*beta.*(tau-gamma).*o;
  o_d = 2*phi.*(delta-D).*o;
  o_tt = 2*beta.*(o+(tau-gamma).*o_t);
  o_td = 2*beta.*(tau-gamma).*o_d;
  o_dd = 2*phi.*(o+(delta-D).*o_d);
  
  i3 = 10:14;
  S3 = sum(s(i3).*o);
  S3_t = sum(s_t(i3).*o+s(i3).*o_t);
  S3_d = sum(s_d(i3).*o+s(i3).*o_d);
  S3_tt = sum(s_tt(i3).*o + 2*s_t(i3).*o_t + s(i3).*o_tt);
  S3_td = sum(s_td(i3).*o + s_t(i3).*o_d + s_d(i3).*o_t + s(i3).*o_td);
  S3_dd = sum(s_dd(i3).*o + 2*s_d(i3).*o_d + s(i3).*o_dd);
  res3 = [S3,S3_t,S3_d,S3_tt,S3_td,S3_dd];

  if max_order > 2
    s_ttt = s_tt.*(t-2)/tau;
    s_ttd = s_tt.*d/delta;
    s_tdd = s_dd.*t/tau;
    s_ddd = s_dd.*(d-2)/delta; 

    S1_ttt = sum(s_ttt(i1));
    S1_ddd = sum(s_ddd(i1)); 
    S1_ttd = sum(s_ttd(i1));
    S1_tdd = sum(s_tdd(i1));
    res1 = [res1,S1_ttt,S1_ttd,S1_tdd,S1_ddd];

    S2_ttt = sum(s_ttt(i2))*e;
    S2_ttd = sum((d(i2)/delta-1).*s_tt(i2))*e;
    S2_tdd = sum(t(i2)/tau.*(s_dd(i2)-2*s_d(i2)+s(i2)))*e;
    S2_ddd = sum(s_ddd(i2)-3*s_dd(i2)+3*s_d(i2)-s(i2))*e;
    res2 = [res2,S2_ttt,S2_ttd,S2_tdd,S2_ddd];

    o_ttt = 2*beta.*(2.*o_t+(tau-gamma).*o_tt);
    o_ttd = 2*beta.*(o_d+(tau-gamma).*o_td);
    o_tdd = 2*phi.*(o_t+(delta-D).*o_td);
    o_ddd = 2*phi.*(2*o_d + (delta-D).*o_dd);

    S3_ttt = sum(s_ttt(i3).*o + 3*s_tt(i3).*o_t + 3*s_t(i3).*o_tt ...
      + s(i3).*o_ttt);
    S3_ttd = sum(s_ttd(i3).*o + s_tt(i3).*o_d + 2*s_td(i3).*o_t ...
      + 2*s_t(i3).*o_td +s_d(i3).*o_tt + s(i3).*o_ttd); 
    S3_tdd = sum(s_tdd(i3).*o + s_dd(i3).*o_t + 2*s_td(i3).*o_d ...
      + 2*s_d(i3).*o_td +s_t(i3).*o_dd + s(i3).*o_tdd);  
    S3_ddd = sum(s_ddd(i3).*o + 3*s_dd(i3).*o_d + 3*s_d(i3).*o_dd ...
      + s(i3).*o_ddd);
  
    res3 = [res3,S3_ttt,S3_ttd,S3_tdd,S3_ddd];
  end
  
  res = res1 + res2 + res3;
end
