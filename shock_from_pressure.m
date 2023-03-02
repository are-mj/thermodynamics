function [u2,T2,v2,M] = shock_from_pressure(th,T1,p1,p2)
  % Calculate shock given pressure on both sides
  %  State 1: undisturbed gas, velocity u1 = 0
  %  State 2: Gas behind shock wave
  % Input:
  %  th: thermodynamic object for the gas. E.g. th = thermo('Air')
  %  T1, p1:  Temperature (K) and pressure (Pa) for State 1
  %  p2: Pressure behind shock (Pa)
  %  z0: Intial vector for newton iteration
  % Output:
  %  u2: Flow velocity behind shock
  %  T2: Temperature behond shock
  %  v2: Molar volume behind shock (m3/kmol)
  %  M:  Shock Mach number
  
  % Reference: Landau, L.D., and E-M. Lifschitz. Fluid Mechanics. 
  %            Pergamon Press, 1959.
%%
  th.Tpcalc(T1,p1);
  % Find ideal gas as starting values for the Newton iteration:
  if nargin < 5
    gamma = th.cp/th.cv;  % Heat cacacity tatio
    R  = th.R;            % Universal gas constant
    Mw = th.Mw;           % Molar mass
    v1_ig = R*T1/p1;         % Molar volume
    c1_ig = sqrt(gamma*p1*v1_ig/Mw);  % Speed of sound
    % Shock Mach number.  L&L equation (85.8): 
    M = sqrt((p2/p1*(gamma+1)+(gamma-1))/(2*gamma)); 
    W_ig = M*c1_ig;
    % Temperature behind shock  L&L (85.9):
    T2_ig = T1*(2*gamma*M^2 - (gamma-1))*((gamma-1)*M^2 + 2) ...
      /(gamma+1)^2/M^2;
    v2_ig = R*T2_ig/p2;
    % Initial values for solving Hugoniot relations:
    z0 = [W_ig;T2_ig;v2_ig];  
  end
  
  % Solve Hugoniot relations  
  
  v1 = th.v;
  c1 = th.c;
  h1 = th.h;
  fun = @(z) residual(th,p1,v1,h1,z,p2);
  z2 = newton(fun,z0);
  W  = z2(1);
  T2 = z2(2);
  v2 = z2(3);
  u2 = (1-v2/v1)*W;
  M = W/c1;  
end

function [f,J] = residual(th,p1,v1,h1,z,p2)
% Residual for Hugoniot realtions
  Mw = th.Mw;
  W  = z(1);
  T2 = z(2);
  v2 = z(3);
  u1_ = -W;
  u1_W = -1;      % du1_/dW
  u2_ = -v2/v1*W;
  u2_v2 = -W/v1;  % du2_/dv2
  u2_W  = -v2/v1; % du2_/dW
  th.Tvcalc(T2,v2);
  h2 = th.h;
  
  f = [0;0;0];
  % Hugoniot conditions:
  % Landau & Lifshitz (82.1), (82.2) and (82.3)
%   f(1) = u2_/v2 - u1_/v1;                       % Conservation of moles
  f(1) = u2_^2/v2 + p2/Mw - u1_^2/v1 - p1/Mw;   % Conservation of momentum
  f(2) = h2 + 0.5*Mw*u2_^2 - h1 - 0.5*Mw*u1_^2; % Conservation of energy
  % Specify p2, calculate u1_
  f(3) = th.p - p2; 
  
  J = zeros(3);
  % Jacobian:
  J(1,1) = 2*u2_/v2*u2_W - 2*u1_/v1*u1_W ;
  J(1,2) = 0;
  J(1,3) = 2*u2_/v2*u2_v2-u2_^2/v2^2;
  
  J(2,1) = Mw*u2_*u2_W - Mw*u1_*u1_W;
  J(2,2) = th.h_T;
  J(2,3) = th.h_v + Mw*u2_*u2_v2;
  
  J(3,2) = th.p_T;
  J(3,3) = th.p_v;
  J(3,1) = 0;
end
