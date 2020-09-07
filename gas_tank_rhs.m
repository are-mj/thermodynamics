function dzdt = gas_tank_rhs(th,z,Ndot,hflow,par)
% tank_rhs(th,T,p,flowfun,par): Time derivatives of T and N for a gas tank
%  To be used for dynamic simulation of filling or emptying
%   th: Thermodynaic object
%   z : [T;N];  T: Gas temperature, N: Molar amount (kmol)
%   w : Mass filling rate (kmol/s)  Ndot < 0 means emptying
%   hflow: Molar (stagnation)enthalpy of gas stream.  
%     If Ndot < 0  hflow is overwritten by h(T,N)
%   par: Tank parameters
% Simplified treatment of temperature profile in the wall:
%    Assumes that par.wallmass takes the same temperature as the  gas
%    Set to some fracrtion of real wall mass.

% Are Mjaavatten, November 2019

  T = z(1);
  N = z(2);
  v = par.V/N;
  Tvcalc(th,T,v);
  Cp_wall = par.wallmass*par.cp_wall;
  Q = par.Area*par.heatcoeff*(T - par.Ta);  % Heat loss
  dNdt = Ndot;
  dTdt = (dNdt*(hflow+v*th.u_v-th.u)-Q)/(N*th.u_T+Cp_wall);
  dzdt = [dTdt;dNdt];
end


  


