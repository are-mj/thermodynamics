    function [ps,ps_T] = psat(T,par)
    % [ps,ps_T] = psat(T,par)  Saturation pressure.  Vectorized
    %  T:    Temperature (K)  Scalar or 1-D array
    %  par:  parameter struct or thermo object
    %  ps:   Saturation pressure (Pa)
    %  ps_T: Derivative (dps/dT)

      if isa(par,'thermo')
          par = par.par;
      end
      T = T(:)';   % Make sure T is a row vector
      theta = 1-T/par.Tc;
      as    = par.as;
      ase   = par.ase;
      if verLessThan('matlab','9.1')
        ps = zeros(1,length(T));
        ps_T = ps;
        for i = 1:length(T)
          xx = par.Tc/T(i).*(as*theta(i).^ase);
          ps(i) = par.pc*exp(xx);
          ps_T(i)  = -ps(i)/T(i)*(xx+as(1) + ...
             as(2:end)*(ase(2:end).*theta(i).^ase(2:end)/theta(i)));
        end       
      else
        xx    = par.Tc./T.*(as*theta.^ase);
        ps    = par.pc*exp(xx);
        ps_T  = -ps./T.*(xx+as(1) + ...
          as(2:end)*(ase(2:end).*theta.^ase(2:end)./theta));
      end
    end