function [T2,v2,W] = compressor(th,T1,p1,p2,eta)
% Outlet conditions from gas compressor
% Input:
%  th:    Thermodynaic object
%  T1:    Inlet temperature (K)
%  p1:    Inlet pressure (Pa)
%  p2:    Outlet pressure
%  eta:   Compressor efficiency (-)
% Output:
%  T2:    Outlet temperature (K)
%  v2-:   Outlet molar volume (m3/kmol)
%  W:     Shaft work (J/kg)

  if ~isa(th,'thermo')
    error('The first argument must be a thermo object')
  end
  if p1 > psat(T1,th)
      error('Initial state must be pure gas, i.e. p1 must be <= psat(T1)')
  end
  th.Tpcalc(T1,p1);  % upstream state
  s1 = th.s;
  h1 = th.h;
  
  % Initialise th with the ideal gas solution:
  gamma = min(th.cp/th.cv,1.6);
  T2_ig = T1*(p2/p1)^((gamma-1)/gamma);
  Tpcalc(th,T2_ig,p2);  
  
  th.pscalc(p2,s1);   % Isentropic compression to pressure p2
  h2 = th.h;          % Molar enthalpy for isentropic process
  h2r = h1 + (h2-h1)/eta;  % Molar enthalpy for real process
  th.phcalc(p2,h2r)
  T2 = th.T;
  v2 = th.v;
  W = (h2r - h1)/th.Mw;  % Shaft work (J(kg)
end