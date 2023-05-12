classdef thermo < handle
  %THERMO: Class for point values of molar thermodynamic variables
  % Calculates molar Helmholtz free energy with partial derivatives with
  % respect to temperature T and molar volume v of orders up to three.
  % This allows explicit expressions for most thermodynamic 
  % variables, a number of which are available as object properties.
  % Type "termo.properties" for a complete list of properties.
  % Type "thermo.methods" for a list of methods.
  %
  % Creator function: th = thermo(species); 
  % Currently available species:  H2, CO2
  % Uses the same H2 model as The NIST Chemistry Webbook
  % (https://webbook.nist.gov/chemistry/fluid/)
  % The results are identical within a factor 5e-5.
  %
  %  The thermo object will contain a consistent set of variables 
  %  when initilalized by either of the following methods
  %   th.Tvcalc(T,v)  % T: temperature (K), v: molar volume (m3/kmol)
  %   th.Tpcalc(T,p)  % T: temperature (K), p: pressure (Pa)
  %
  % Please note that the native variables are T and v, so Tvcalc requires
  % no iterations.  Tpcalc involves iterative searches for T and v, 
  % and will converge in most cases.  
  %
  % Are Mjaavatten, January 2020
  % Mjaavatten Consulting, Norway.  are@mjaavatten.com 

  % Lowercase letters are used for all intensive properties except
  % temperature
  
  properties
    species;   % Name of current species (e.g. 'H2')    
    R;         % Universal gas constant (J/kmol/K)
    Mw;        % Molar mass (kg/kmol)
    Tc;        % Critical temperature (K)
    pc;        % critical pressure
    rhoc;      % critical density (kmol/m3)
    vc;        % Critical molar volume (m3/kmol)
    par;       % Model parameters
    Tt;        % Triple point temperature (K)
    pt;        % Triple point pressure (Pa)
    max_order; % Highest level of partial derivatives 
               %  Default: 2, Maximum: 3
  end  

  properties  % Thermodynamic variables
    helmholtz % Function handke tor thermodynmic model
    T;        % Temperature (K)
    v;        % Volume (m3/kmol)
    f;        % Helmholtz free energy (J/kmol)
    f_T;      % df/dT = -s
    f_v;      % df/dv = -p
    f_TT;     % d2^f/dT^2 = -ds/dT
    f_Tv;     % d2^f/dTdv = d2^f/dvdT = -ds/dv = -dp/dT
    f_vv;     % d^2f/dv^2  = -dp/dv
    f_TTT;    % d^3f/dT^3
    f_TTv;    % d^3f/dT^2dv
    f_Tvv;    % d^3f/dTdv^2
    f_vvv;    % d^3/dv*3;
    p;        % Pressure (Pa)
    p_T;      % dp/dT
    p_v;      % dp/dv;
    p_TT;     % d2p/dT^2
    p_Tv;     % d2p/dTdv
    p_vv;     % d2p/dv^2
    u;        % Internal energy (j/kmol)
    u_T;      % du/dT
    u_v;      % du/dv
    cv;       % Heat capacity at constant volume (J/kmol/K)
    cp;       % Heat capacity at constant pressure (J/kmol/K)
    s;        % Entropy (J/kmol/K)
    s_T;      % ds/dT
    s_v;      % ds/dv
    h;        % Enthalpy (J/kmol)
    h_T;      % dh/dT
    h_v;      % dh/dv
    g;        % Gibbs energy (J/kmol)
    c;        % Speed of sound (m/s)
    c_T;      % dc/dT
    c_v;      % dc/dv
    jt;       % Joule-Thompson coefficent (dT/dp at constant h)    
  end

  methods
    function th = thermo(species,Mw,gamma)
    % THERMO: Create new thermodynamic object
    % Usage: th = thermo(species);
    %     or th = thermo(species,Mw,gamma)  (For species 'ig')
      if nargin < 1
        error(['You must specify the chemical species, ', ...
          'e.g. th = thermo(''H2'')']);
      end
      species_list = {'H2','orthoH2','paraH2','CO2','H2O','N2','O2' ...
        ,'Ar','Air','ig'};
      sp = find(strcmpi(species_list,species),1);
      if isempty(sp)
        fprintf('''%s'' not found. Currently available species are:\n  '...
          ,species);
        for i = 1:length(species_list)-1
          fprintf(' ''%s'', ',species_list{i})
        end
        fprintf('and ''%s''\n ',species_list{end});
        error('Species not found')
      end
      if strcmp(species,'ig')
        if nargin <2
          Mw = input('Ideal gas molar mass (kg/kmol)? ');
          gamma = input('Ideal gas heat capacity ratio? ');
        end
        th.par = parameters_ig(Mw,gamma);
      else
        th.par = feval(['parameters_',species]);
      end
        % Unpack general parameters:
      th.species = th.par.species;
      th.R    = th.par.R;      % Universal gas constant (J/(kmol K)
      th.Tc   = th.par.Tc;     % Critical temperature (K)
      th.pc   = th.par.pc;     % Critical pressure (Pa)
      th.rhoc = th.par.rhoc;   % Critical molar density (kmol/m3)
      th.Mw   = th.par.Mw;     % Molar mass (kg/kmol)
      th.vc   = th.par.vc;     % Critical molar volume (m3/kmol)
      th.Tt   = th.par.Tt;     % Triple point temperature (K)
      th.pt   = th.par.pt;     % Triple point temperature (K)
      th.max_order = 2;
      th.helmholtz = @helmholtz;
    end
        
    function Tvcalc(th,T,v)
    % TVCALC: Initialise thermodynamic object from temperature and volume
    %  Explicit calculation, no iterations
    % Input:
    %   T :   Temperature (K)
    %   v :   Molar volume (m3/kmol)
    
      th.T = T;
      th.v = v;
      res = th.helmholtz(T,v,th.par,th.max_order);
      th.f = res(1);
      th.f_T = res(2);
      th.f_v = res(3);
      th.f_TT = res(4);
      th.f_Tv = res(5);
      th.f_vv = res(6);

      % derived properties:
      th.p = -th.f_v;       % Pressure
      th.p_T = -th.f_Tv;
      th.p_v = -th.f_vv;
      th.s = -th.f_T;       % Entropy
      th.s_T = -th.f_TT;
      th.s_v = -th.f_Tv;
      
      th.u = th.f + T*th.s;   % Internal energy
      th.u_T = th.f_T + th.s + T*th.s_T;
      th.u_v = th.f_v +T*th.s_v;
      th.h = th.u + th.p*v;   % Enthalpy
      th.h_T = th.u_T + th.p_T*v;
      th.h_v = th.u_v + th.p_v*v + th.p;
      th.g   = th.f + th.p*th.v;    % Molar Gibb'senergy

      th.cv = th.u_T;
      th.cp = th.cv + T*th.f_Tv^2/th.f_vv;
      th.c = sqrt(v^2/th.Mw*(th.f_vv-th.f_Tv^2/th.f_TT)); %speed of sound
      th.jt = -(T*th.f_Tv/th.f_vv+v)/th.cp;  % Joule-Thompson coeff.

      if th.max_order < 3  
        if ~isempty(th.f_TTT)
          % Erase third-order derivatives to avoid errors later
          th.f_TTT = [];
          th.f_TTv = [];
          th.f_Tvv = [];
          th.f_vvv = [];
          th.p_TT  = [];
          th.p_Tv  = [];
          th.p_vv  = [];
          th.c_T   = [];
          th.c_v   = [];          
        end
      else
        th.f_TTT = res(7);
        th.f_TTv = res(8);
        th.f_Tvv = res(9);
        th.f_vvv = res(10);
        th.p_TT = -th.f_TTv;
        th.p_Tv = -th.f_Tvv;
        th.p_vv = -th.f_vvv; 	
        % Derivatives of c:
        th.c_T = v^2/th.Mw*(th.f_Tvv-2*th.f_Tv/th.f_TT*th.f_TTv ...
          + th.f_Tv^2/th.f_TT^2*th.f_TTT)/2/th.c;
    	  th.c_v = (2*v/th.Mw*(th.f_vv-th.f_Tv^2/th.f_TT) ...
          + v^2/th.Mw*(th.f_vvv-2*th.f_Tv/th.f_TT*th.f_Tvv ...
          + th.f_Tv^2/th.f_TT^2*th.f_TTv))/2/th.c;          
      end
    end

    function Tpcalc(th,T,p,v0)
    % Tpcalc(T,p,v0) : Solves p(T,v) = p by Newton's method
	  % v0: optional starting point for iteration. Default: Ideal gas value
      if nargin < 4
        v0 = th.R*T/p;  % Ideal gas
      end
      MaxIter = 25;
      vv = v0;
      for i = 1:MaxIter
        th.Tvcalc(T,vv)
        ff = th.p-p;
        if abs(ff) < 0.01
          return
        end
        J = th.p_v;
        dv = ff/J;
        while vv - dv < 0
          dv = dv/2;
        end
        vv = vv - dv;
      end
      error('No convergence for T = %6.2f, p = %10.0f',T,p);
    end

    function phcalc(th,p,h)
    % phcalc(p,h): Initialise thermodynamic object from pressure and enhalpy
    %  Solves (p(T,v) = p, h(T,v) = h by Newton's method
    % Input:
    %   p:   Pressure (Pa)
    %   h:   Molar enthalpy (J/kmol) 
    % Initial T and v are calculated from the current thermodynamic object
    % assuming constant Joule-Thompson coefficient.
      if isempty(th.T)
        error(['Please initialise the thermo object state ' ,...
          'using either Tvcalc or Pvcalc']);
      end    
      MaxIter = 25;
      Tol = 1e-10;
      % Estimate new T from current jt:
      T1 = th.T + th.jt*(p-th.p);
      th.Tpcalc(T1,p);
      x = [T1;th.v];
      dx = 1000;
      for i = 1:MaxIter
        if norm(dx)<Tol
            break
        end
        th.Tvcalc(x(1),x(2));
        fun = [th.p-p;th.h-h];
        J = [th.p_T th.p_v;th.h_T th.h_v];
        dx = J\fun;
        x = x-dx;
      end  
      if i >= MaxIter
        error('No convergence for T = %6.2f, p = %10.0f',...
            x(1),x(2));
      end
    end      
        
    function pscalc(th,p,s)
    % pscalc(p,s): Initialise thermodynamic object from pressure 
    % and molar entropy
    % Solves p(T,v) = p; s(T,v) = s
    % Input:
    %   p:   Pressure (Pa)
    %   s:   Molar entropy (J/kmol(K) 
    % Initial T and v is taken from the thermodynamic object
    % Use Tvcalc ot Tpcalc to define a good starting point    
      if isempty(th.T)
        error(['Please initialise the thermo object state ' ,...
          'using either Tvcalc or Pvcalc']);
      end        
      MaxIter = 25;
      Tol = 1e-10;
      x = [th.T;th.v];
      dx = 1000;      
      for i = 1:MaxIter
        if norm(dx)<Tol
            break
        end
        th.Tvcalc(x(1),x(2));
        fun = [th.p-p;th.s-s];
        J = [th.p_T th.p_v;th.s_T th.s_v];
        dx = J\fun;
        x = x-dx;
      end     
        if i >= MaxIter
            error('No convergence for T = %6.2f, p = %10.0f',...
                x(1),x(2));
        end
    end

    function [ps,vl,vv,ps_T] = saturation(th,T)
    % saturation(T): Values along the vapour/liquid saturation line
    %   T:    Temperature (K)
    %  Output:
    %   ps:   Saturation pressure (Pa)
    %   vl:   Saturated liquid molar volume (M3/kmol)
    %   vv:   Saturated vapour molar volume (M3/kmol)
    %   ps_T: Derivative of pd (dps/dT)
      ps = NaN;
      vl = NaN;
      vv = NaN;
      ps_T = NaN;  
      if ~isfield(th.par,'as')
        error('Saturation data for %s are missing',th.species)
      end
      if T > th.Tc
        warning('Temperature T = %0.2f is above the critical temperature tC = %.2f',T,th.Tc)
        return
      end
      theta  = 1-T/th.Tc;
      as = th.par.as;
      ase = th.par.ase;
      xx    = th.Tc/T*as*theta.^ase;
      ps    = th.pc*exp(xx);
      ps_T  = -ps/T*(xx+as(1) + ...
        as(2:end)*(ase(2:end).*theta.^ase(2:end)/theta));
      vl   = th.vc/(th.par.bs*theta.^th.par.bse);
      vv   = th.vc*exp(th.par.cs*theta.^th.par.cse);  
    end

    function s = thermprops(th)
      s.T = th.T;
      s.v = th.v;
      s.f = th.f;
      s.f_T = th.f_T;
      s.f_v = th.f_v;
      s.f_TT = th.f_TT;
      s.f_Tv = th.f_Tv;
      s.f_vv = th.f_vv;
      s.f_TTT = th.f_TTT;
      s.f_TTv = th.f_TTv;
      s.f_Tvv = th.f_Tvv;
      s.f_vvv = th.f_vvv;
      s.p = th.p;
      s.p_T = th.p_T;
      s.p_v = th.p_v;
      s.p_TT = th.p_TT;
      s.p_Tv = th.p_Tv;
      s.p_vv = th.p_vv;
      s.u = th.u;
      s.u_T = th.u_T;
      s.u_v = th.u_v;
      s.cv = th.cv;
      s.cp = th.cp;
      s.s = th.s;
      s.s_T = th.s_T;
      s.s_v = th.s_v;
      s.h = th.h;
      s.h_T = th.h_T;
      s.h_v = th.h_v;
      s.g   = th.g;
      s.c = th.c;
      s.jt = th.jt;
    end
  end
  
  methods(Static)

    function properties
    % PROPERTIES: Displays a list of all propetries of the thermo class
      names = {'R','species','Mw','Tc','pc','rhoc','Tt','par','max_order'};
      text = {
        'Universal gas constant (J/kmol/K)'	
        'Name of current species (e.g., ''H2'')'
        'Molar mass (kg/kmol)'
        'Critical temperature (K)'
        'critical pressure (Pa)'
        'critical density (kmol/m3)'
        'Triple point temperature, K'
        'Parameters from Leachman et al. (2005)'
        'Maximum derivative order (default: 2)'
          };
      fprintf('%12s\n','Parameters:')
      fprintf('%10s   %-50s\n','Name','Description');
      for i = 1:length(names)
          fprintf('%10s : %-50s\n',names{i},text{i});
      end  

      names = {'T','v','f','f_T','f_v','f_TT','f_Tv','f_vv','p','p_T',...
        'p_v','p_TT','p_Tv','p_vv','u','u_T','u_v','cv','cp','s','s_T',...
        's_v','h','h_T','h_v','g','c','jt'};
      text = {
        'Temperature (K)'
        'Volume (m3/kmol)'
        'Helmholtz free energy (J/kmol)'
        'df/dT = -s'
        'df/dv = -p'
        'd2^f/dT^2 = -ds/dT'
        'd2^f/dTdv = d2^f/dvdT = -ds/dv = -dp/dT'
        'd^2f/dv^2  = -dp/dv'
        'Pressure (Pa)'
        'dp/dT'
        'dp/dv;'
        'd2p/dT^2'
        'd2p/dTdv'
        'd2p/dv^2'
        'Internal energy (J/kmol)'
        'du/dT'
        'dh/dv'
        'Heat capacity at constant volume (J/kmol/K)'
        'Heat capacity at constant pressure (J/kmol/K)'
        'Entropy (J/kmol/K)'
        'ds/dT'
        'ds/dv'
        'Enthalpy (J/kmol)'
        'dh/dT'
        'dh/dv'
        'Gibbs free energy (J/kmol)'
        'Speed of sound (m/s)'
        'Joule-Thompson coefficent (dT/dp at constant h)'
        };
      names3 = {'f_TTT','f_TTv','f_Tvv','f_vvv','c_T','c_v'};      
      text3 = {
        'd^3f/dT^3'
        'd^3f/dT^2dv'
        'd^3f/dTdv^2'
        'd^3/dv*3'
        'dc/T'
        'dc/dv'
        };
      fprintf('\n Properties:\n')
      fprintf('%10s   %-50s\n','Name','Description');
      for i = 1:length(names)
          fprintf('%10s : %-50s\n',names{i},text{i});
      end
      fprintf('\n\tProperties that are calculated if max_order = 3\n')
      fprintf('%10s   %-50s\n','Name','Description');
      for i = 1:length(names3)
          fprintf('%10s : %-50s\n',names3{i},text3{i});
      end      
    end

    function methods
    % METHODS: Displays a list of methods of the thermo class
      names = {'thermo','Tvcalc','Tpcalc','phcalc',...
        'pscalc','saturation','thermprops','properties','methods'};
      text = {'Creator function for thermodynamic object'
        'Calculate thermo, given temperature and molar volume'
        'Calculate thermo, given temperature and pressure'
        'Calculate thermo, given pressure and molar enthalpy'
        'Calculate thermo, given pressure and molar entropy'
        'Saturation pressure, saturated liquid and vapour thermo objects'
        'Return struct with all thermodynamin properties'
        'List thermo properties'
        'List thermo methods'
          };  
      fprintf('\n Methods:\n')
      fprintf('%10s   %-50s\n','Name','Description');
      for i = 1:length(names)
        fprintf('%10s : %-50s\n',names{i},text{i});
      end 
      fprintf('Note:\n')
      fprintf(['\tphcalc and pscalc are iterative metods that ', ...
        'require iniial values for T and v.\n']);
      fprintf(['\tthese are taken from the thermo object. Therefore ', ...
        'the object must first \n'])
      fprintf(['\tbe initalised, e.g. using Tvcalc', ...
        ' or Tpcalc.\n'])
      fprintf('\nType: help thermo.<method> for more help\n') 
    end
  end    
end

