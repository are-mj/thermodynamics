function [T,ps,vv,vl,dh] = vapour_pressure(species,T)
% Pressure and enthalpy of evaporation along the saturation curve
%  species:  'H2','CO2' of 'H2O'
%  T:   Temperature array (K). If missing, T = linspace(th.Tt,th.Tc)
%  ps:  Vapour pressure (Pa)
%  vl:  molar volume of saturated liquid  (m3/kmol)
%  vv:  molar volume of saturated vapour  (M3/kmol)
%  dh:  enthalpy of vaporisation (J/kmol)

  plotting = true;
  th = thermo(species);
  if nargin < 2
    T = linspace(th.Tt,th.Tc);
  end
  if length(T) < 4
    plotting = false;
  end
  ps = zeros(size(T));
  vl = zeros(size(T));
  vv = zeros(size(T));
  dh = zeros(size(T));
  for i = 1:length(T)
    [ps(i),vl(i),vv(i)] = th.saturation(T(i));
    th.Tvcalc(T(i),vl(i));
    hl = th.h;
    th.Tvcalc(T(i),vv(i));
    hv = th.h;    
    dh(i) = hv-hl;
  end
  
  if plotting
    figure;
    subplot(211);
    plot(T,ps*1e-5);
    ylabel('p (bar)');
    title(sprintf('Saturation pressure of %s',species));
    grid on
    subplot(212);
    plot(T,dh*1e-3/th.Mw);
    xlabel('T (K)');
    ylabel('\Deltah (kJ/kg)');
    title('Enthalpy of evaporation')
    grid on
  end
    
    