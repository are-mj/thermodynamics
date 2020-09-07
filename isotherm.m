function [p,vl,vv,ps] = isotherm(th,T,v)
% [p,vl,vv,ps] = isotherm(th,T,v): Isotherm at temperature T.
%   T: Temperature (K)
%   v: Array of molar volumes (m3/kmol)
% Output:
%   p:  Array of pressures at constant temperature T, as function of v
%  If T < Tc:
%   vl,vv: Liquid and vapour molar volumes at ps.
%   ps: Saturation pressure at T
%
% Example:
%   th = thermo('H2');
%   v = logspace(log10(0.0275),1,1000); T = 25;
%   [p,vl,vv,ps] = isotherm(th,T,v);
%   figure; semilogx(v,p*1e-5,[vl,vv],ps*1e-5*[1,1],'k',[vl,vv],ps*1e-5*[1,1],'*r')
%   ylabel bar; xlabel m^3/kmol
%   legend('Model','Saturation pressure','Saturation volumes')
%   axis([0,10,-10,10])
%   grid on
%   title('Isotherm for H_2 at 25K')

  v = v(:)';
  npoints = length(v);
  p = zeros(1,npoints);
  for i = 1:npoints
    th.Tvcalc(T,v(i));
    p(i) = th.p;
  end
  vl = NaN;
  vv = NaN;
  ps = NaN;
  if T < th.Tc
    [ps,vl,vv] = th.saturation(T);
  end
end
  
