function [t,Z,pA,pB] = gas_tank_system(th,T_up,p_up,T_down,p_down,Ta,w)
% Simulate system with constant flowrate between two gas tanks
%   T_up,p_up:  Initial temperature and pressure of upstream reservoir
%   T_down,p_down: Initial temperature and pressure of downstream reservoir
%   Ta       :  Ambient temperature (K)
%   w        : (Constant) Mass flow rate
%
% Example:
% Ta = 273.15+15;
% T_up = Ta;
% p_up = 350e5;
% T_down = Ta;
% p_down = 20e5;
% w = 0.2;  % kg/s
% gas_tank_system(T_up,p_up,T_down,p_down,Ta,w);
  
  % Upstream tank parameters
  tankA.wallmass = 500;    % kg
  tankA.cp_wall = 100;     % J/kg/K 
  tankA.V = 5;             % m3
  tankA.Area = 2;          % m2
  tankA.heatcoeff  = 0;    % W/m2/K  (Adiabatic tank)
  Tpcalc(th,T_up,p_up);
  tankA.N0 = tankA.V/th.v; % Initial kmoles
  tankA.Ta = Ta;           % Ambient temperature (K)
  hfill = th.h;
  
  % Downstream tank 
  tankB.wallmass = 500;    % kg
  tankB.cp_wall = 100;     % J/kg/K 
  tankB.V = 5;             % m3
  tankB.Area = 2;          % m2
  tankB.heatcoeff  = 0;    % W/m2/K
  Tpcalc(th,T_down,p_down);
  tankB.N0 = tankB.V/th.v; % Initial kmoles
  tankB.Ta = Ta;           % Ambient temperature (K)
  
  Z0 = [T_up;tankA.N0;T_down;tankB.N0];
  
  Ndot = w/th.Mw;
  fun = @(t,Z) tanksystem_rhs(Z,th,Ndot,hfill,tankA,tankB);
  stopfun = @(t,Z) no_dp(Z,th,tankA,tankB);
  opt = odeset('events',stopfun);
  [t,Z] = ode45(fun,[0,1000],Z0,opt);
  
  TA = Z(:,1);
  NA = Z(:,2);
  TB = Z(:,3);
  NB = Z(:,4);
  
  nn = length(t);
  pA = zeros(nn,1);
  pB = zeros(nn,1);
  for i = 1:length(t)
    Tvcalc(th,TA(i),tankA.V/NA(i));
    pA(i) = th.p;
    Tvcalc(th,TB(i),tankB.V/NB(i));
    pB(i) = th.p;    
  end
end

function dZdt = tanksystem_rhs(Z,th,Ndot,hfill,tankA,tankB)
  dZdt(1:2,1) = gas_tank_rhs(th,Z(1:2),-Ndot,hfill,tankA);
  dZdt(3:4,1) = gas_tank_rhs(th,Z(3:4),Ndot,hfill,tankB);
end

function [value,isterminal,direction] = no_dp(Z,th,tankA,tankB)
% Terminates integration when the pressure difference between tanks is zero
  Tvcalc(th,Z(1),tankA.V/Z(2));
  pA = th.p;
  Tvcalc(th,Z(3),tankB.V/Z(4));
  pB = th.p;
  value = pA-pB; 
  isterminal = 1;
  direction = 0;
end
