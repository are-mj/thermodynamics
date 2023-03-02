function twophase_region(species)
%  Plots the envelope of the two-phase region in the v-p plane

switch species
  case {'H2','paraH2','orthoH2'}
    ax = [0.02,20,0,20];
    v = logspace(log10(2e-2),log(20),500);
    Th = 35;  % Isotherm temperatures
    Tl = 28;
  case 'CO2'
    ax = [0.03,4,0,100];
    v = logspace(log10(2e-2),log(20),500);
    Th = 310;  % Isotherm temperatures
    Tl = 295;   
  case 'H2O'
    ax = [0.018,20,0,300];
    v = logspace(log10(0.018),log(20),500);
    Th = 655;  % Isotherm temperatures
    Tl = 620;          
  otherwise
    error('Currently available species are: H2, paraH2, orthoH2, CO2, H2O')
end
    
th = thermo(species);

% Calculate saturation pressure and liquid and vapour molar volumes
% for a list of temperatures.  Plot saturation envelope (ps vs. v)
n = 500;T = linspace(th.Tt,th.Tc,n)';
ps = zeros(n,1);v_liq = zeros(n,1);v_vap = zeros(n,1);
for i = 1:n-1
  [ps(i),v_liq(i),v_vap(i)] = th.saturation(T(i));
end
ps(n) = th.pc;
ps = ps*1e-5;
v_liq(n) = th.vc;
v_vap(n) = th.vc;
figure;
semilogx([v_liq;flipud(v_vap)],[ps;flipud(ps)],'b');
title(sprintf('Two-phase region and isotherms for %s',species))
ylabel('Pressure (bar)')
xlabel('Volume (m^3/kmol)')

% Generate isotherms for Tl, Tc and Th
p_l = isotherm(th,Tl,v)*1e-5;     % Subclritical
pc = isotherm(th,th.Tc,v)*1e-5;   % Critical
p_h = isotherm(th,Th,v)*1e-5;     % Supercritical

% Plot isotherm for Tl.  This is unphysocal in the two-phase region.
hold on;
plot(v,p_l,':k');
% Find the saturated pressure and molar volumess at Tl:
[pl,vl,vv] = th.saturation(Tl);
% Remove unphysical values and add values for vl and vv:
liq_index = find(v<vl);
vap_index = find(v>vv);
real_v = [v(liq_index),vl,vv,v(vap_index)];
real_p_l = [p_l(liq_index),[pl,pl]*1e-5,p_l(vap_index)];
plot(real_v,real_p_l,'k');  % PLot real isotherm
plot([vl,vv],pl*1e-5*[1,1],'k','marker','*'); % Mark saturation points

% Critical isotherm
plot(v,pc,'r');
plot(1/th.rhoc,th.pc*1e-5,'*r');

% Supercritical isotherm
plot(v,p_h,'--k');
lines = flipud(get(gca,'children'));
subset = lines([1,2,4,5,6,7]);
leg2 = sprintf('%dK model pressure',Tl);
leg3 = sprintf('%dK real isotherm',Tl);
leg6 = sprintf('%dK model pressure',Th);
legend(subset,{'Two-phase envelope',leg2,leg3,'Critical isotherm',...
  'Critical point',leg6})
axis(ax)
end

