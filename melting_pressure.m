function pm = melting_pressure(T,par)
% MELTING_PRESSURE: Melting pressure as function of T

beta = T/par.Tt-1;
pm = par.pt*(1+par.melta(1)*beta + par.melta(2)*beta.^2);

end

