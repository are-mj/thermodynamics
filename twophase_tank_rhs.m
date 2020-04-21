function dzdt = twophase_tank_rhs(th,t,z,par)
% Time derivatives when gas or liquid flows into or out of a twophase tank
%   th:  thermo object
%   z:  state vector:
%      N:  Total tank content (kmol)
%      T:  Temperature (K)
%      vl: Liquid molar volume at saturation (m3/kmol)
%      vv: Vapour molar volume at saturation (m3/kmol)
%   par: parameters struct:
%      Vtank   : Volume of tank (m3)
%      flowfun : User-supplied function giving the molar flowrate and  
%                molar enthalpy of the flow entering or leaving the tank
%                 [Ndot,h] = flowfun(th,t,z,par)
%                 Ndot mut be < 0 for flow out of the tank
%
% The true state variables are U (or T) and N.  To avoid a complex 
%  Differential-Algebraic Equation system we use N,T,vl and vv as 
%  state variables.

  [Ndot,hflow] = par.flowfun(th,t,z,par);
  
  N  = z(1); % Tank molar content (kmol)
  T  = z(2); % Temperature (K)
  vl = z(3); % Liquid molar volume at equilibrium
  vv = z(4); % Vapour molar volume at equilibrium
  
  [~,~,~,ps_T] = th.saturation(T);
  
  % Liquid phase
  th.Tvcalc(T,vl)
  ul = th.u;
  hl = th.h;
  
  dvldT = (ps_T-th.p_T)/th.p_v;
  duldT = th.u_T + th.u_v*dvldT;

  % Vapour phase
  th.Tvcalc(T,vv)
  uv = th.u; 
  
  dvvdT = (ps_T-th.p_T)/th.p_v;
  duvdT = th.u_T + th.u_v*dvvdT;
  
  dNdt = Ndot;
  dUdt = Ndot*hflow + par.Q;
  
  Nl = (N*vv-par.Vtank)/(vv-vl);
  Nv = N - Nl;
  dNldN = vv/(vv-vl);
  dNvdN = 1-dNldN;
  dNldT = N*dvvdT/(vv-vl)-(N*vv-par.Vtank)/(vv-vl)^2*(dvvdT-dvldT);
  dNvdT = -dNldT;
  
  dUdT = dNldT*ul + Nl*duldT + dNvdT*uv + Nv*duvdT;
  dUdN = dNldN*ul+dNvdN*uv;
  dTdt = (dUdt-dUdN*dNdt)/dUdT;
  
  dzdt = [dNdt;dTdt;dvldT*dTdt;dvvdT*dTdt];
end
