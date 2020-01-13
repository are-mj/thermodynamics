function [t,T,N,p] = fill_gas_tank(th,tank,T0,p0,T_up,p_up,Ta,w,p_final)
% Simulate the filling of a gas tank from a constant gas reservoir
%  Input:
%   th       : Thermodynamic object
%   tank     : Paramteter struct for the tank
%   T0,p0    : Initial temperature (K) and pressure (Pa) in tank
%   T_up,p_up:  Stagnation temperature and pressure of upstream reservoir
%   Ta       :  Ambient temperature (K)
%   w        : (Constant) Mass flow rate
%   p_final  :  Stop simulation when the tank pressure reaches p_final
% Output:
%   t        : time
%   T        : Tank temperature (K)
%   N        : Tank contents (kmol)
%   p        : Tank pressure (Pa)

% January 2020, Are Mjaavatten

  th.Tpcalc(T0,p0);
  tank.N0 = tank.V/th.v; % Initial kmoles
  tank.Ta = Ta;           % Ambient temperature (K)
  
  Tpcalc(th,T_up,p_up);
  hfill = th.h;    % Stagnation enthalpy of feed
  
  Tpcalc(th,T0,p0)
  v0 = th.v;
  N0 = tank.V/v0;   % Number of moles in tank at start
  z0 = [T0;N0];
  
  Ndot = w/th.Mw;
  fun = @(t,z) gas_tank_rhs(th,z,Ndot,hfill,tank);
  stopfun = @(t,z) tank_full(z,th,tank,p_final);
  opt = odeset('events',stopfun);
  [t,z] = ode45(fun,[0,1000],z0,opt);
  
  T = z(:,1);
  N = z(:,2);
  
  nn = length(t);
  p = zeros(nn,1);
  for i = 1:length(t)
    Tvcalc(th,T(i),tank.V/N(i));
    p(i) = th.p;
  end
end

function [value,isterminal,direction] = tank_full(z,th,tank,p_final)
% Terminates integration when the pressure dfference between tanks is low
  Tvcalc(th,z(1),tank.V/z(2));
  value = p_final-th.p; 
  isterminal = 1;
  direction = 0;
end
