% CO2twophase_region.m
% Script for creating a volume-pressure plot of the two-phase region
clear
th = thermo('CO2');

% 
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
title('Two-phase region and isotherms for CO_2')
ylabel('Pressure (bar)')
xlabel('Volume (m^3/kmol)')
ax = [0.03,3,0,100];

v = logspace(log10(2e-2),log(20),2000);
p295 = isotherm(th,295,v)*1e-5;
pc = isotherm(th,th.Tc,v)*1e-5;
p310 = isotherm(th,310,v)*1e-5;
hold on;
plot(v,p295,':k')

% The real isotherm has constant pressure in the two-phase region
% Intersection points isotherm/two-phase envelope:
[vil,pil] = polyxpoly(v_liq,ps,v,p295);
vil = min(vil);
[viv,piv] = polyxpoly(v_vap,ps,v,p295);
viv = max(viv);
piv = min(piv);

il = find(v < vil,1,'last');  % Highest v in liquid region
ih = find(v > viv,1,'first'); % Lowest v in vapour region
v_295 = [v(1:il),vil,viv,v(ih:end)];       
p_295 = [p295(1:il),piv,piv,p295(ih:end)];  % Real isotherm

plot(v_295,p_295,'k')
plot(v,pc,'r')
plot(th.vc,th.pc*1e-5,'*r')
plot(v,p310,'--k')
legend('Two-phase envelope','295K model pressure','295K real isotherm',...
  'Critical isotherm','Critical point','310K isotherm')

axis(ax)



