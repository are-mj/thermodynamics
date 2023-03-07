function demo_tank_leak
% Simulation of leak from the vapour phase in a tank containing liquid H2
% Demonstrates calculation of vapour/liquid equilibrium as well as
% choked and subsonic flow.
%
% NOTE:  Liquid will condense in the escaping gas.  This effect is ignored.
%  Corresponds to the "frozen equilibrium" assumption for two-phase orifice 
%  flow.  Real flow rates may differ from the ones calculated here!

  % Parameters:
  par.Vtank = 1;             % 1 m3
  par.A_leak = pi*0.005^2;   % 10 mm leak
  par.Q = 0;                 % Adiabatic tank
  par.p_a = 101325;          % Ambient pressure (Pa)
  
  % Initial conditions
  T0 = 28;   % K
  N0 = 10;   % kmol
  
  % Simulation parameters
  tmax = 400;  % nax simulation time (s)
  Tol = 1e-4;  % Relative tolerance in ode45

  % Flow parameters
 
  th = thermo('H2');  % thermodynamic object

  [~,vl,vv] = th.saturation(T0);
  z0 = [N0;T0;vl;vv];
  
  z1 = z0';  % to ensure correct starting point in flow is not choked
  t1 = 0;

  % Check if flow is initially choked:
  warning('off','CHOKED_FLOW:RANGE')
  choked = sonic_gap(th,0,z0,par)>0;
  if choked  % Stage 1: Sonic flow
    % Integrate as long as flow is sonic:
    stopfun = @(t,z) sonic_gap(th,t,z,par);
    par.flowfun = @choked_leak;
    opt = odeset('events',stopfun,'RelTol',Tol);
    dzdt = @(t,z) twophase_tank_rhs(th,t,z,par);
    [t1,z1] = ode23(dzdt,[0,tmax],z0,opt);
  end
  % Continue with subsonic flow as long as the pressure difference is
  % positive:
  par.flowfun = @subsonic_leak;
  stopfun = @(t,z) pressure_gap(th,t,z,par);
  opt = odeset('events',stopfun,'RelTol',Tol);
  dzdt = @(t,z) twophase_tank_rhs(th,t,z,par);
  [t2,z2] = ode23(dzdt,[t1(end),tmax],z1(end,:)',opt);

  % combine and extract variables
  t = [t1;t2];
  z = [z1;z2];
%   U  = z(:,1);
  N  = z(:,1);
  T  = z(:,2);
  vl = z(:,3);
  vv = z(:,4);
  x = (vv-par.Vtank./N)./(vv-vl);  % Liquid fraction

  % Calculate pressures in liquid and vapour phases.
  % Start out as equal, but accumulating numerical errors may cause values
  % to drift (Numerical drift)
  pl = zeros(size(t));
  for i = 1:length(t)
    th.Tvcalc(T(i),vl(i));
    pl(i) = th.p;
  end
  pv = zeros(size(t));
  for i = 1:length(t)
    th.Tvcalc(T(i),vv(i));
    pv(i) = th.p;
  end

  % Calculate conditions at exit:
  c1eak1 = zeros(size(t1));
  pleak1 = zeros(size(t1));
  Tleak1 = zeros(size(t1));
  for i = 1:length(t1)
    [~,~,Tleak1(i),pleak1(i),c1eak1(i)] = choked_leak(th,t,z1(i,:)',par);
  end
  cleak2 = zeros(size(t2));
  uleak2 = zeros(size(t2));
  Tleak2 = zeros(size(t2));
  for i = 1:length(t2)
    [~,~,Tleak2(i),uleak2(i),cleak2(i)] = subsonic_leak(th,t,z2(i,:)',par);
  end
  Tleak = [Tleak1;Tleak2];
  cleak = [c1eak1;cleak2];
  uleak = [c1eak1;uleak2];
  pleak = [pleak1;par.p_a*ones(size(t2))];
  
  warning('on','CHOKED_FLOW:RANGE')
  
  figure;
  subplot(221);
  plot(t,[T,Tleak]);
  title('Temperatures');legend('Tank','Exit')
  ylabel K
  subplot(222); 
  plot(t,th.Mw*[N,N.*x,N.*(1-x)]);
  legend('Total','Liquid','Vapour','location','se');
  title('Mass of H_2 in tank'); ylabel kg
  subplot(223); plot(t,[pv,pleak]*1e-5); title Pressure
  xlabel('Time (s)');ylabel bar; legend('Tank','Exit')

  subplot(224); plot(t,uleak,t,cleak); xlabel('Time (s)')
  title('Velocities at exit')
  ylabel('m/s');legend('Flow velocity','Speed of sound','location','sw')
  xlabel('Time (s)');

  fprintf('Total leak:  %10.3f kg\n',(N0-N(end))*th.Mw)
end

function [value,isterminal,direction] = sonic_gap(th,t,z,par)
  % The difference between the sonic exit pressure and the
  % ambient flow.  The flow is choked as long as this value is positive
  % Event function for ODE solver
  [~,~,~,pstar] = choked_leak(th,t,z,par);
  value = pstar - par.p_a;
  isterminal = 1;
  direction = -1;
end 

function [Ndot,hflow,Tstar,pstar,cstar] = choked_leak(th,t,z,par)
% Flow variables at choking flow from vapour phase of tank
%  z: State vector for tank
%  Ndot:  Molar flow rate (kmol/s)
%  hflow: Molar stagnstion enthalpy in flow
%  Tstar: Exit temperature (K)
%  pstar: Exit pressure (Pa)
%  cstar: Exit speed of sound (and flow speed)
  Tvap = z(2);
  vvap  = z(4);
  th.Tvcalc(Tvap,vvap)
  pvap = th.p;
  hflow = th.h;  
  [Tstar,pstar,cstar,vstar] = choked_flow(th,Tvap,pvap);
  Ndot = -par.A_leak*cstar/vstar;
end

function [Ndot,h0,T_out,u_out,c_out] = subsonic_leak(th,~,z,par)
% Flow rate for leak from the vapour phase
  MaxIter = 25;
  Tolx = 1e-5;
  T0 = z(2);
  v0 = z(4);  % Vapour molar volume
  th.Tvcalc(T0,v0)
  h0 = th.h;
  s0 = th.s;
  % Subsonic flow:
  x = [T0;v0];
  dx = inf;
  for Iter = 1:MaxIter
    if norm(dx) < Tolx
      break
    end    
    T = x(1);
    v = x(2);
    th.Tvcalc(T,v);
    f = [th.s - s0
         th.p - par.p_a];
    J = [th.s_T,th.s_v
         th.p_T,th.p_v];   
    dx = -J\f;
    x = x + dx;
  end 
  if Iter >= MaxIter
    error('No convergence in vapour_leak')
  end
  v_out = th.v;  % Molar volune at outlet
  u_out = sqrt(2/th.Mw*max(h0-th.h,0));  % Outlat flow velocity
  c_out = th.c;
  T_out = th.T;
  Ndot = -par.A_leak*u_out/v_out;
end
  
function [value,isterminal,direction] = pressure_gap(th,~,z,par)
  % The pressure difference between the tank and the ambient
  % Event function for ODE solver
  th.Tvcalc(z(2),z(4));
  deltap = 10; % Pa
  value = th.p - (par.p_a +deltap);
  isterminal = 1;
  direction = -1;
end

