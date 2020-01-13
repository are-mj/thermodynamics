% twophase_region.m
clear
th = thermo('H2');
n = 200;T = linspace(th.Tt,th.Tc,n)';
ps = zeros(n,1);v_liq = zeros(n,1);v_vap = zeros(n,1);
for i = 1:n-1
  [ps(i),liq,vap] = th.saturation(T(i));
  v_liq(i) = liq.v;
  v_vap(i) = vap.v;
end
ps(n) = th.pc;
ps = ps*1e-5;
v_liq(n) = th.vc;
v_vap(n) = th.vc;

figure;semilogx([v_liq;flipud(v_vap)],[ps;flipud(ps)],'b')
title('Two-phase region and isotherms for H_2')
ylabel('Pressure (bar)')
xlabel('Volume (m^3/kmol)')
ax = [0.02,20,0,20];

v = logspace(log10(2e-2),log(20),2000);
p28 = isotherm(th,28,v)*1e-5;
pc = isotherm(th,th.Tc,v)*1e-5;
p35 = isotherm(th,35,v)*1e-5;
hold on;
plot(v,p28,':k')

% The real isotherm has constant pressure in the two-phase region
% Intersection points isotherm/two-phase envelope:
[vil,pil] = polyxpoly(v_liq,ps,v,p28);
vil = min(vil);
[viv,piv] = polyxpoly(v_vap,ps,v,p28);
viv = max(viv);
piv = min(piv);
il = find(v < vil,1,'last');  % Highest v in liquid region
ih = find(v > viv,1,'first'); % Lowest v in vapour region
v_28 = [v(1:il),vil,viv,v(ih:end)];       
p_28 = [p28(1:il),piv,piv,p28(ih:end)];  % Real isotherm

plot(v_28,p_28,'k')
plot(v,pc,'r')
plot(1/th.rhoc,th.pc*1e-5,'*r')
plot(v,p35,'--k')
legend('Two-phase envelope','28K model pressure','28K real isotherm',...
  'Critical isotherm','Critical point','35K isotherm')

axis(ax)


