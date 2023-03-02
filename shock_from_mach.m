function [p2,T2,v2,u2] = shock_from_mach(th,T1,p1,M)
  % Calculate shock given shock Mach number
  %  State 1: undisturbed gas, velocity u1 = 0
  %  State 2: Gas behind shock wave
  % Input:
  %  th: thermodynamic object for the gas. E.g. th = thermo('Air')
  %  T1, p1:  Temperature (K) and pressure (Pa) for State 1
  %  M:  Shock Mach number
  % Output:
  %  p2: Pressure behind shock (Pa)  
  %  T2: Temperature behond shock
  %  v2: Molar volume behind shock (m3/kmol)
  %  u2: Flow velocity behind shock
  
  % Reference: Landau, L.D., and E-M. Lifschitz. Fluid Mechanics. 
  %            Pergamon Press, 1959.

  th.Tpcalc(T1,p1);
  c1 = th.c;
  W = c1*M;
  %% Find ideal gas as starting values for the Newton iteration:
  gamma = th.cp/th.cv;  % Heat cacacity ratio
  R  = th.R;            % Universal gas constant
  v1_ig = R*T1/p1;      % Molar volume
  % L&L (85.9):
  T2_ig = T1*(2*gamma*M^2 - (gamma-1))*((gamma-1)*M^2 + 2)/(gamma+1)^2/M^2;
  v2_ig = v1_ig*((gamma-1)*M^2+2)/(gamma+1)/M^2;
  z0 = [T2_ig;v2_ig];
  
  %% Solve Hugoniot relations  
  v1 = th.v;
  h1 = th.h;
  fun = @(z) residual(th,p1,v1,h1,W,z);
  z2 = newton(fun,z0);
  T2 = z2(1);
  v2 = z2(2);
  u2 = W*(1 - v2/v1);
  p2 = th.p;
end

function [f,J] = residual(th,p1,v1,h1,W,z)
% Residual for Hugoniot realtions
  Mw = th.Mw;
  T2 = z(1);
  v2 = z(2);
  u1_ = -W;
  u2_ = -v2/v1*W;
  th.Tvcalc(T2,v2);
  p2 = th.p;
  h2 = th.h;
  
  u2_v2 = u1_/v1; % du2_/dv2

  f = [0;0];
  % Hugoniot conditions:
  % Landau & Lifshitz (82.1), (82.2) and (82.3)
%   f(1) = u2_/v2 - u1_/v1;                       % Conservation of moles
  f(1) = u2_^2/v2 + p2/Mw - u1_^2/v1 - p1/Mw;   % Conservation of momentum
  f(2) = h2 + 0.5*Mw*u2_^2 - h1 - 0.5*Mw*u1_^2; % Conservation of energy
  
  % Jacobian:
  J = zeros(2);
  J(1,1) = th.p_T/Mw;
  J(1,2) = 2*u2_/v2*u2_v2-u2_^2/v2^2 + th.p_v/Mw;
  J(2,1) = th.cp;
  J(2,2) = th.h_v + Mw*u2_*u2_v2;
end
