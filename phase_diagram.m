function phase_diagram(species)
  % Draws the bundaries between solid, liquid and gas phases in T-p diagram
  
  switch species
    case 'H2O'
      par = parameters_H2O;
      melt = @(T) melting_pressure(T,par);
      sublim = @(T) H2Osublimation_pressure(T,par);
      Tm = [linspace(260,272),linspace(272,par.Tt,500)];
      Tsub = linspace(220,par.Tt);
      ax = [220,650,2e-5,1.5e3];
    case 'CO2'
      par = parameters_CO2;
      melt = @(T) melting_pressure(T,par);
      sublim = @(T) CO2sublimation_pressure(T,par);
      Tm = linspace(par.Tt,250);
      Tsub = linspace(160,par.Tt); 
      ax = [160,310,0.03,1500];
    otherwise
      fprintf('Phase diagram implemented only for H2O and CO2\n');
      return
  end
  Ts = linspace(par.Tt,par.Tc,500);
  ps = psat(Ts,par);
  pm = melt(Tm);
  psub = sublim(Tsub);
  figure; semilogy(Ts,ps*1e-5,Tm,pm*1e-5,Tsub,psub*1e-5);
  xlabel('T (K)');
  ylabel('p (bar)')
  legend('Liquid-gas','Solid-liquid','Solid-gas','location','southeast')
  title(sprintf('Phase diagram for %s',species))
  grid on
  axis(ax);
end

function psub = CO2sublimation_pressure(T,par)
% Sublimation pressure as function of temperature
% CO2 only
  % Eq. (3.12) in Span-Wagner
  theta = 1-T/par.Tt;
  a = par.sublima;
  e = par.sublime;  
  if verLessThan('matlab','9.1')
    xx = zeros(1,length(T));
    for i = 1:length(T)
      xx(i) = par.Tt/T(i)*(a*theta(i).^e);
    end       
  else
    xx    = par.Tt./T.*(a*theta.^e);
  end
  psub  = par.pt*exp(xx);
end

function psub = H2Osublimation_pressure(T,par)
% Sublimation pressure as function of temperature
% Wagner, W. et al., J. Phys. Chem. Ref. Data 40.4 (2011)
  theta = T/par.Tt; 
  a = par.sublima;
  b = par.sublime;
  xx = 0;
  for i = 1:numel(a)
    xx = xx + a(i)*theta.^b(i);
  end
  psub = par.pt*exp(xx./theta);
end