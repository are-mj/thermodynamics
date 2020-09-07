function phase_diagram(species)
  % Draws the bundaries between solid, liquid and gas phases in T-p diagram
  switch species
    case 'H2O'
      par = H2Oparameters;
      sublim = @(T) H2Osublimation_pressure(T,par);
      melt = @(T) H2Omelting_pressure(T,par);
      Tm = [linspace(260,272),linspace(272,par.Tt,500)];
      Tsub = linspace(220,par.Tt);
      ax = [220,650,2e-5,1.5e3];
    case 'CO2'
      par = CO2parameters;
      sublim = @(T) CO2sublimation_pressure(T,par);
      melt = @(T) CO2melting_pressure(T,par);  
      Tm = linspace(par.Tt,250);
      Tsub = linspace(160,par.Tt); 
      ax = [160,310,0.03,1500];
    otherwise
      fprintf('Phase diagram not implemented for %s\n',species);
      return
  end
  Ts = linspace(par.Tt,par.Tc);
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

function pm = CO2melting_pressure(T,par)
% CO2MELTING_PRESSURE: Melting pressure as function of T for CO2
beta = T/par.Tt-1;
pm = par.pt*(1+par.melta(1)*beta + par.melta(2)*beta.^2);
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

function pm = H2Omelting_pressure(T,par)
% H2OMELTING_PRESSURE: Melting pressure as function of T for H2O
  theta = T/par.Tt;
  a = par.melta;
  e = par.melte;
  pm = par.pt*(1+a(1)*(1-theta.^e(1)) + a(2)*(1-theta.^e(2)));
end

function psub = H2Osublimation_pressure(T,par)
% Sublimation pressure as function of temperature
% CO2 only
  % Eq. (3.12) in Span-Wagner
  theta = T/par.Tt;
  a = par.sublima;
  e = par.sublime;  
  xx = a(1)*(1-theta.^e(1))+a(2)*(1-theta.^e(2));
  psub = par.pt*exp(xx);
