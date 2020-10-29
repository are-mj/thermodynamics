function [T2,p2,u2,T3,T5,p5] = shocktube(T1,p1,species1,T4,p4,species4)
% Calculate shock tube with rigorous thermodynamics
% Input:
%    T1, p1:  Temperature (K) and pressure (Pa) for undisturbed driven gas
%    species1:  Thermo object or species name for driven gas. E.g. 'Air'
%    T4, p4:  Temperature (K) and pressure (Pa) for undisturbed driver gas
%    species4:  Thermo object or species name for driver gas. E.g. 'H2'
% Output:
%  T2,p2,u2:  Temperature, molar volume and velocity in region 2 
%             (driven gas behind shock). 
%  T3:        Temperature in region 3 (expanded driver gas)
%             Note that u3 = u2, p3 = p2
%  T5,p5:     Temperature and pressure in shock reflected from end of 
%             driven gas section
% Shock tube zones:
%    1: Undisturbed driven gas 
%    2: Shocked, driven gas
%    3: Driver gas between rarefaction fan and contact surface
%    4: Undisturbed driver gas
%    5: Reflected shock
%  
% Usage example:
%   [T2,v2,u2,T3,v3] = shocktube(300,1e5,'Air',300,100e5,'H2')
%   th2 = thermo('Air');th2.Tvcalc(T2,v2);p2 = th2.p;

% Are Mjaavatten, September 2020

% Solution strategy:
%  Calculate pressure p3 and velocity u3 for an epansion wave from T4,p4 
%   to p3 by solving du/dv = c/v for constant entropy and inital conditions
%   v0 = v4, u(v4) = 0, T(0) = T4,  (s = s4 = constant)
%   Stop integration when p3 = shock_from_velocity(th1,T1,p1,u1,u2);
%  Finally calculate the reflected shoch

  if isa(species1,'thermo')
    % Useful for debugging using perfect gas
    th1 = species1;
    th4 = species4;
  else
    th1 = thermo(species1);
    th4 = thermo(species4);
  end
  th1.Tpcalc(T1,p1);
  th4.Tpcalc(T4,p4);
  % Initial values
  z0 = [0;th4.T];   
  v0 = th4.v;
  fun = @(v,z) Riemann(th4,v,z);
  opt = odeset('events',@(v,z) stopfun(th1,th4,T1,p1,v,z));
  [v,z] = ode45(fun,v0*[1,10],z0,opt);
  T3 = z(end,2);
  v3 = v(end); 
  th4.Tvcalc(T3,v3);
  u2 = z(end,1);    % Common flow velocity in zones 2 and 3
  u1 = 0;
  [p2,T2] = shock_from_velocity(th1,T1,p1,u1,u2);
  % Reflected shock:
  u5 = 0;
  [p5,T5] = shock_from_velocity(th1,T2,p2,u2,u5);
end

function dzdv = Riemann(th,v,z)
  % Find the Riemann invariant: intergral(c(v)/v*dv)-u
  % by solving du/dv = c/v for constant entropy
  T = z(2);
  th.Tvcalc(T,v);
  dudv = th.c/v;
  dTdv = -th.s_v/th.s_T;
  dzdv = [dudv;dTdv];
end

function [value, isterminal, direction] = stopfun(th1,th4,T1,p1,v3,z3)
  u1 = 0;   % Undisturbed gas in front of shock
  u2 = z3(1);
  T3 = z3(2);
  th4.Tvcalc(T3,v3)
  p3 = th4.p;
  p2 = shock_from_velocity(th1,T1,p1,u1,u2);
  value      = p2-p3;
  isterminal = 1;   % Stop the integration
  direction  = 0;
end 


   
  
