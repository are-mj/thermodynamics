function [p2,T2,v2,M] = shock_from_velocity(th,T1,p1,u1,u2)
  %  Calculate shock given velocities on both sides
  %  th: thermodynamic object
  %  T1: Temperature in front of shock (K)
  %  p1: Pressure in front of shock (Pa)
  %  u1: Flow velocity in front of shock (m/s)
  %  u2: Flow velocity behind shock (m/s)
  %
  %  p2,T2,v2: Pressure, temperature and molar volume behind shock
  %  M : Shock Mach number
  
  if abs(u2-u1) < 1e-5  % No shock
    p2 = p1;
    T2 = T1;
    M = 0;
    v2 = th.v;
    return
  end
    
  th.Tpcalc(T1,p1);
  v1 = th.v;
  c1 = th.c;
  T1 = th.T;
  % Initial guess (Ideal gas values):
  gamma = th.cp/th.cv;
  uu = u2*(gamma+1)/4/c1;
  M10 = uu + sqrt(uu^2+1);  % Mach number for perfect gas
  W0 = M10*c1*sign(u2-u1);  % Shock velocity for perfect gas
  T20 = T1 + 0.5*th.Mw*((u1-W0)^2-(u2-W0)^2)/th.cp;  % Energy conservation
  v20 = v1*(u2-W0)/(u1-W0); % Mass conservation
  z0 = [W0;T20;v20];        % Initial values for Newton solver
  fun = @(z) residual(th,u1,u2,T1,v1,z);
  z = newton(fun,z0);
  W = z(1);
  T2 = z(2);
  v2 = z(3);
  th.Tvcalc(T2,v2);
  p2 = th.p;
  M = W/c1;
end

function [f,J] = residual(th,u1,u2,T1,v1,z)
  th.Tvcalc(T1,v1);
  p1 = th.p;
  h1 = th.h;
  W  = z(1); 
  T2 = z(2);
  v2 = z(3);
  u1_ = u1-W;
  u2_ = u2-W;
  th.Tvcalc(T2,v2);
  p2 = th.p;
  h2 = th.h;
  Mw = th.Mw;
  f = [0;0;0];
  f(1) = u2_/v2 - u1_/v1;                       % Conservation of moles
  f(2) = u2_^2/v2 + p2/Mw - u1_^2/v1 - p1/Mw;   % Conservation of momentum
  f(3) = h2 + 0.5*Mw*u2_^2 - h1 - 0.5*Mw*u1_^2; % Conservation of energy
  % Jacobian J(i,j) = df(i)/dz(j) :
  J = zeros(3);
  J(1,1) = 1/v1-1/v2;
  J(1,3) = -u2_/v2^2;
  J(2,1) = -2*(u1_/v1-u2_/v2);
  J(2,2) = th.p_T/Mw;
  J(2,3) = -u2_^2/v2^2 + th.p_v/Mw;
  J(3,1) = Mw*(u1_-u2_);
  J(3,2) = th.h_T;
  J(3,3) = th.h_v;
end