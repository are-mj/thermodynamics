function pm = melting_pressure(T,par)
% MELTING_PRESSURE: Melting pressure as function of T

theta = T/par.Tt;
if strcmp(par.species,'H2O')
  xx = 0;
  for i = 1:numel(par.melta)
    xx = xx + par.melta(i)*(1-theta.^par.melte(i));
  end
  pm = par.pt*(1+xx);
% Ref: Wagner, W. et al., J. Phys. Chem. Ref. Data 40.4 (2011)
else
  beta = theta-1;
  pm = par.pt*(1+par.melta(1)*beta.^par.melte(1) ...
    + par.melta(2)*beta.^par.melte(2));
% Ref: Giordano, V. et al., J.Chem. Phys, 
end