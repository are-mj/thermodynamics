function res_out = helmholtz(T,v,par,max_order)
% Helmholtz molar free energy (a) and partial derivatives  
% Input:
%   T:   Temperature (K)
%   v:   Molar volume (m3/kmol)
%   par: Thermodynamic model parameter struct
%   max_order: Highest order of partial derivatives. Default: 2, Max: 3
% Output:
%   res_out = [a,a_T,a_v,a_TT,a_Tv,a_vv]  (a_xy = d^2a/dxdy)
%    if max_order = 3: res_out = [res_out,a_TTT,a_TTv,a_Tvv,a_vvv]

% March 2020, Are Mjaavatten

  if nargin < 4
    max_order = 2;
  end
  
  % Reduced temperature and density:
  tau = par.Tc/T;
  delta = par.vc/v;
  tau_T = -tau/T;
  delta_v = -delta/v;
  delta_vv = -2*delta_v/v;
  
  if isfield(par,'phi_ig')
    res_ig = par.phi_ig(tau,delta,par,max_order);
  else
    res_ig = phi_ig(tau,delta,par,max_order);
  end
  if isfield(par,'n')
    res_r = phi_r(tau,delta,par,max_order);
    res = res_ig+res_r;
  else  % ideal gas
    res = res_ig;
  end

  % Reduced Helmholts free energy and deriventive w.r.t tau (t) and delta (d)
  phi = res(1);
  phi_t = res(2);
  phi_d = res(3);
  phi_tt = res(4);
  phi_td = res(5);
  phi_dd = res(6);
  
  R = par.R;
  a = R*T*phi;
  a_T = R*(phi + T*phi_t*tau_T);
  a_v = R*T*phi_d*delta_v;
  % a_TT = R*(2*phi_t*tau_T + T*phi_tt*tau_T^2+T*phi_t*tau_TT);
  %  But terms 1 and 3 cancel, so:
  a_TT = R*T*phi_tt*tau_T^2;
  a_Tv = R*(phi_d*delta_v + T*phi_td*tau_T*delta_v);
  a_vv = R*T*(phi_dd*delta_v^2 + phi_d*delta_vv);
  res_out = [a,a_T,a_v,a_TT,a_Tv,a_vv];
  if max_order > 2
    tau_TT = -2*tau_T/T;
    delta_vvv = -3*delta_vv/v;
    phi_ttt = res(7);
    phi_ttd = res(8);
    phi_tdd = res(9);    
    phi_ddd = res(10);
    a_TTT = R*(phi_tt*tau_T^2 + T*phi_ttt*tau_T^3 ...
      + 2*T*phi_tt*tau_T*tau_TT);
    a_TTv = R*T*phi_ttd*tau_T^2*delta_v;
    a_Tvv = R*(phi_dd*delta_v^2 + phi_d*delta_vv ...
      + T*phi_tdd*tau_T*delta_v^2 + T*phi_td*tau_T*delta_vv);
    a_vvv = R*T*(phi_ddd*delta_v^3 + 3*phi_dd*delta_v*delta_vv ...
      + phi_d*delta_vvv);
    res_out = [res_out,a_TTT,a_TTv,a_Tvv,a_vvv];
  end
end

function res = phi_ig(tau,delta,par,max_order)
  a = par.ig_a;
  if isfield(par,'ig_b')
    N = numel(par.ig_a);
    theta = par.ig_b;
    A = a(4:N)';
    B = theta';
    y = exp(-tau*B)';
  else
    A     = 0;
    B     = 0;
    y     = 0;
  end
  phi = log(delta) + a(1) + a(2)*tau + a(3)*log(tau) + A*log(1-y);
  phi_t = a(2) + a(3)/tau + A.*B*(y./(1-y));
  phi_d = 1/delta;
  phi_tt = -a(3)/tau^2 - (A.*B.^2)*(y./(1-y).^2);
  phi_dd = -delta^(-2);
  phi_td = 0;
  res = [phi,phi_t,phi_d,phi_tt,phi_td,phi_dd];
  if max_order > 2
    phi_ttt = 2*a(3)/tau^3+(A.*B.^3)*(y.*(1+y)./(1-y).^3);
    phi_ddd = 2*delta^(-3);
    res = [res,phi_ttt,0,0,phi_ddd];
  end
end

function res = phi_r(tau,delta,par,max_order)
% Reduced residual Helmholtz free energy a with derivatives
% Create phi_r as a sum of contributions phi_r = S1 + S2 + S3 + S4
  pos = par.sections;
  dpos = diff(pos);
  sections = length(pos);
  if sections >2
    n = par.n(1:pos(3));
  else
    n = par.n(1:pos(2));
  end
  t = par.t;
  d = par.d;
  c = par.c;
   
  s = n.*delta.^d.*tau.^t;
  s_t = s.*t/tau;
  s_d = s.*d/delta;
  s_tt = s_t.*(t-1)/tau;
  s_td = s_t.*d/delta;
  s_dd = s_d.*(d-1)/delta; 

  i1 = 1:pos(1);
  S1 = sum(s(i1));
  S1_t = sum(s_t(i1));
  S1_d = sum(s_d(i1));  
  S1_tt = sum(s_tt(i1));
  S1_dd = sum(s_dd(i1));
  S1_td = sum(s_td(i1));
  res1 = [S1,S1_t,S1_d,S1_tt,S1_td,S1_dd];
  
  i2 = (pos(1)+1):pos(2);
  e = exp(-delta.^c);
  e_d = -c.*delta.^(c-1).*e;
  e_dd = e_d.*((c-1)/delta+e_d./e);
  S2 = s(i2)'*e;
  S2_t = s_t(i2)'*e;
  S2_d = s_d(i2)'*e+s(i2)'*e_d;
  S2_tt = s_tt(i2)'*e;
  S2_td = s_td(i2)'*e + s_t(i2)'*e_d;
  S2_dd = s_dd(i2)'*e + 2*s_d(i2)'*e_d + s(i2)'*e_dd;
  res2 = [S2,S2_t,S2_d,S2_tt,S2_td,S2_dd];
  
  if sections > 2
    alpha = par.alpha;
    beta3 = par.beta(1:dpos(2));
    epsilon = par.epsilon;
    gamma = par.gamma;     
    
    o = exp(-alpha.*(delta-epsilon).^2 - beta3.*(tau-gamma).^2);
    o_t = -2*beta3.*(tau-gamma).*o;
    o_d = -2*alpha.*(delta-epsilon).*o;
    o_tt = -2*beta3.*(o+(tau-gamma).*o_t);
    o_td = -2*beta3.*(tau-gamma).*o_d;
    o_dd = -2*alpha.*(o+(delta-epsilon).*o_d);

    i3 = (pos(2)+1):pos(3);
    S3 = sum(s(i3).*o);
    S3_t = sum(s_t(i3).*o+s(i3).*o_t);
    S3_d = sum(s_d(i3).*o+s(i3).*o_d);
    S3_tt = sum(s_tt(i3).*o + 2*s_t(i3).*o_t + s(i3).*o_tt);
    S3_td = sum(s_td(i3).*o + s_t(i3).*o_d + s_d(i3).*o_t + s(i3).*o_td);
    S3_dd = sum(s_dd(i3).*o + 2*s_d(i3).*o_d + s(i3).*o_dd);
    res3 = [S3,S3_t,S3_d,S3_tt,S3_td,S3_dd];
  else
    res3 = zeros(1,6);
  end

  if sections > 3
    a = par.a;
    b = par.b;
%     beta3 = par.beta(1:dpos(2));
    beta4 = par.beta(dpos(2)+(1:dpos(3)));  
    A = par.A(1);  % exploit that A(1:3) = constant
    B = par.B;
    C = par.C;
    D = par.D; 
    
    if abs(delta-1)<1e-15
      resGamma = ones(dpos(3),1)*[0,0,0,Inf,0,0];
      if max_order > 2
        resGamma = [resGamma,ones(dpos(3),1)*[Inf,Inf,0,0]];
      end
    else
      x = ((delta-1)^2).^(1/2./beta4);
      x_d = 1./beta4*(delta-1).*abs(delta-1).^(1./beta4-1);
      x_dd = x_d.*(1./beta4-1)/(delta-1);

      theta = 1 - tau + A*x;
      theta_d = A*x_d;
      theta_dd = A.*x_dd;

      y = ((delta-1).^2).^a;
      y_d = 2*a.*y/(delta-1);
      y_dd = (2*a-1).*y_d/(delta-1);

      Delta = theta.^2 + B.*y;
      Delta_t = -2*theta;
      Delta_d = 2*theta.*theta_d + B.*y_d;
      Delta_tt = 2*ones(dpos(3),1);
      Delta_td = -2*theta_d;
      Delta_dd = 2*theta_d.^2+2*theta.*theta_dd+B.*y_dd; 

      Phi = Delta.^b;
      Phi_t = b.*Delta.^(b-1).*Delta_t;
      Phi_d = b.*Delta.^(b-1).*Delta_d;
      bb = b.*(b-1).*Delta.^(b-2);
      Phi_tt = bb.*Delta_t.^2 +  b.*Delta.^(b-1).*Delta_tt;
      Phi_td = bb.*Delta_t.*Delta_d + b.*Delta.^(b-1).*Delta_td;
      Phi_dd = bb.*Delta_d.^2 +  b.*Delta.^(b-1).*Delta_dd;  

      Psi = exp(-C*(delta-1)^2-D*(tau-1)^2);
      Psi_t = -2*D*(tau-1).*Psi;
      Psi_d = -2*C*(delta-1).*Psi;
      Psi_tt = -2*D.*(Psi+(tau-1)*Psi_t);
      Psi_td = -2*D*(tau-1).*Psi_d;
      Psi_dd = -2*C.*(Psi+(delta-1)*Psi_d);

      Gamma = Phi.*Psi*delta;
      Gamma_t = (Phi_t.*Psi + Phi.*Psi_t)*delta;
      Gamma_d = (Phi_d.*Psi + Phi.*Psi_d)*delta + Phi.*Psi;
      Gamma_tt =(Phi_tt.*Psi+2*Psi_t.*Phi_t+Phi.*Psi_tt)*delta;
      Gamma_td = (Phi_td.*Psi+Phi_t.*Psi_d+Phi_d.*Psi_t+Phi.*Psi_td)*delta ...
         + Gamma_t/delta;
      Gamma_dd = (Phi_dd.*Psi+2*Psi_d.*Phi_d+Phi.*Psi_dd)*delta ...
        + 2*(Phi_d.*Psi + Phi.*Psi_d);

      resGamma = [Gamma,Gamma_t,Gamma_d,Gamma_tt,Gamma_td,Gamma_dd];
    end
    n = par.n(pos(3)+1:pos(4))';
    res4 = n*resGamma;
  else
    res4 = zeros(1,6);
  end
  
  if max_order > 2
    s_ttt = s_tt.*(t-2)/tau;
    s_ttd = s_tt.*d/delta;
    s_tdd = s_td.*(d-1)/delta;
    s_ddd = s_dd.*(d-2)/delta;     
    
    i1 = 1:7; 
    S1_ttt = sum(s_ttt(i1));
    S1_ttd = sum(s_ttd(i1));
    S1_tdd = sum(s_tdd(i1));
    S1_ddd = sum(s_ddd(i1));
    res1 = [res1,S1_ttt,S1_ttd,S1_tdd,S1_ddd];
    e_ddd = e_dd.*((c-1)/delta+e_d./e)-e_d.*(c-1).*(1+c.*delta.^c)/delta^2;

    S2_ttt = s_ttt(i2)'*e;
    S2_ttd = s_ttd(i2)'*e + s_tt(i2)'*e_d;
    S2_tdd = s_tdd(i2)'*e + 2*s_td(i2)'*e_d + s_t(i2)'*e_dd;
    S2_ddd = s_ddd(i2)'*e + 3*s_dd(i2)'*e_d + 3*s_d(i2)'*e_dd ...
           + s(i2)'*e_ddd;
    res2 = [res2,S2_ttt,S2_ttd,S2_tdd,S2_ddd];
    
    if sections > 2
      o_ttt = -2*beta3.*(2*o_t+(tau-gamma).*o_tt);
      o_ttd = -2*beta3.*(o_d+(tau-gamma).*o_td);
      o_tdd = -2*alpha.*(o_t+(delta-epsilon).*o_td);
      o_ddd = -2*alpha.*(2*o_d+(delta-epsilon).*o_dd);

      S3_ttt = sum(s_ttt(i3).*o + 3*s_tt(i3).*o_t + 3*s_t(i3).*o_tt ...
             + s(i3).*o_ttt);
      S3_ttd = sum(s_ttd(i3).*o +s_tt(i3).*o_d ...
                   + 2*s_td(i3).*o_t + 2*s_t(i3).*o_td ...
                   + s_d(i3).*o_tt + s(i3).*o_ttd);
      S3_tdd = sum(s_tdd(i3).*o + s_dd(i3).*o_t ...
                   + 2*s_td(i3).*o_d + 2*s_d(i3).*o_td...
                   + s_t(i3).*o_dd + + s(i3).*o_tdd);
      S3_ddd = sum(s_ddd(i3).*o + 3*s_dd(i3).*o_d + 3*s_d(i3).*o_dd ...
             + s(i3).*o_ddd);   
      res3 = [res3,S3_ttt,S3_ttd,S3_tdd,S3_ddd];
    else
      res3 = zeros(1,10);
    end
    if sections> 3
      if abs(delta-1)>1e-15
        x_ddd = (1/beta4-1)/(delta-1)*(x_dd-x_d/(delta-1));
        y_ddd = (2*a-1)/(delta-1).*(y_dd-y_d/(delta-1));
        theta_ddd = A*x_ddd;

        Delta_ttt = zeros(dpos(3),1);
        Delta_ttd = zeros(dpos(3),1);
        Delta_tdd = -2*theta_dd;
        Delta_ddd = 4*theta_d.*theta_dd + 2*theta_d.*theta_dd ...
          + 2*theta.*theta_ddd + B.*y_ddd; 
        Phi_ttt = b.*(b-1).*(b-2).*Delta.^(b-3).*Delta_t.^3 ...
                + 3*bb.*Delta_t.*Delta_tt + b.*Delta.^(b-1).*Delta_ttt;          
        Phi_ttd = b.*(b-1).*(b-2).*Delta.^(b-3).*Delta_t.^2.*Delta_d ...
                + 2*bb.*Delta_t.*Delta_td + bb.*Delta_tt.*Delta_d ...
                + b.*Delta.^(b-1).*Delta_ttd;
        Phi_tdd = b.*(b-1).*(b-2).*Delta.^(b-3).*Delta_t.*Delta_d.^2 ...
                + 2*bb.*Delta_d.*Delta_td + bb.*Delta_t.*Delta_dd ...
                + b.*Delta.^(b-1).*Delta_tdd;   
        Phi_ddd = b.*(b-1).*(b-2).*Delta.^(b-3).*Delta_d.^3 ...
                + 3*bb.*Delta_d.*Delta_dd + b.*Delta.^(b-1).*Delta_ddd;

        Psi_ttt = -2*D.*(2*Psi_t + (tau-1)*Psi_tt);
        Psi_ttd = -2*D.*(Psi_d + (tau-1)*Psi_td);
        Psi_tdd = -2*C.*(Psi_t + (delta-1)*Psi_td);
        Psi_ddd = -2*C.*(2*Psi_d + (delta-1)*Psi_dd);   

        Gamma_ttt = (Phi_ttt.*Psi+2*Phi_tt.*Psi_t+Phi_t.*Psi_tt ...
          +Phi_tt.*Psi_t+2*Phi_t.*Psi_tt+Phi.*Psi_ttt)*delta;
        Gamma_ttd = (Phi_ttd.*Psi+2*Phi_td.*Psi_t+Phi_d.*Psi_tt ...
          + Phi_tt.*Psi_d+2*Phi_t.*Psi_td+Phi.*Psi_ttd)*delta ...
          + Phi_tt.*Psi+2*Phi_t.*Psi_t+Phi.*Psi_tt;
        Gamma_tdd = (Phi_tdd.*Psi+2*Psi_td.*Phi_d+Phi_t.*Psi_dd ...
          + Phi_dd.*Psi_t+2*Psi_d.*Phi_td+Phi.*Psi_tdd)*delta ...
          + 2*(Phi_td.*Psi + Phi_t.*Psi_d + Phi_d.*Psi_t + Phi.*Psi_td);
        Gamma_ddd = (Phi_ddd.*Psi+2*Psi_dd.*Phi_d+Phi_d.*Psi_dd ...
          + Phi_dd.*Psi_d+2*Psi_d.*Phi_dd+Phi.*Psi_ddd)*delta ...
          + Phi_dd.*Psi+2*Psi_d.*Phi_d+Phi.*Psi_dd ...
          + 2*(Phi_dd.*Psi + Phi_d.*Psi_d + Phi_d.*Psi_d + Phi.*Psi_dd);

        resGamma = [resGamma,Gamma_ttt,Gamma_ttd,Gamma_tdd,Gamma_ddd];
      end
      res4 = n*resGamma;
    else
      res4 = zeros(1,10);
    end
  end
  res = res1 + res2 + res3 +res4;
end