function [T,p,c,v] = choked_flow(th,T0,p0)
% State of choked gas flow, given upstream stagnation conditions.
%  Assumes isentropic process and no condensation
% Input:  
%   th   thermo object
%   T0   Stagnation temperature (K)
%   p0   Stagnation pressure (Pa)
% Output:  Variables at choked condiditons
%   T    Temperature    (K)
%   p    Pressure (Pa)
%   c    Speed of sound (m/s)
%   v    Molar volume  (m3/kmol)

% January 2020, Are Mjaavatten
    
    ps = th.saturation(T0);
%     if p0 > ps
%       warning('CHOKED_FLOW:RANGE',...
%         'Calculation not robust in the liquid region')
%     end
    h0 = th.h;
    s0 = th.s;    
    % Initialize with ideal gas solution:
%     gamma = th.cp/th.cv;
    % The ideal gas formula is very wrong near the saturation line
    % This modifiaction seems to work fine for H2:
%     gamma = min(gamma,1.6);  
%     Tstar = T0*2/(gamma+1);
%     vstar = th.v*(T0/Tstar)^(1/(gamma+1));
%     th.Tvcalc(Tstar,vstar);
%     cstar = th.c;
    th.max_order = 3;  % Need third order derivatives for the Jacobian
    fun = @(x) hscfun(th,x,h0,s0);
%     x0 = [Tstar;vstar;cstar];
    x0=[T0;th.v;th.c];
    x = newton(fun,x0,[zeros(3,1),Inf(3,1)]);
    th.max_order = 2;
    % If this simple Newtonsolver fails, try fsolve (Optimization tlbx):
    % x = fsolve(fun,x0,optimoptions('fsolve','display','off'));
    T = x(1);
    v = x(2);
    c = x(3);
    p = th.p;
end

function [fun,J] = hscfun(th,x,h1,s1)
% Residuals of stagnation enthalpy, entropy and flow speed - speed of sound
% hscfun(x,th,h1,s1) = [0;0;0] at choking
    T = x(1);
    v = x(2);
    u = x(3);
    Mw = th.Mw;
    th.Tvcalc(T,v);
    fun = zeros(3,1);
    fun(1) = th.h + 0.5*Mw*u^2 - h1;  % Stagnation enthalpy conserved
    fun(2) = th.s - s1;               % Isentropic flow
    fun(3) = th.c^2 - u^2;            % Flow at speed of sound
    dhdT = -(v*th.f_Tv + T*th.f_TT);
    dhdv = -(v*th.f_vv + T*th.f_Tv);
    f_vv = th.f_vv;
    f_Tv = th.f_Tv;
    f_TT = th.f_TT;
    dc2dT = v^2/Mw*(th.f_Tvv-2*f_Tv/f_TT*th.f_TTv + f_Tv^2/f_TT^2*th.f_TTT);
    dc2dv = 2*v/Mw*(f_vv-f_Tv^2/f_TT) ...
        + v^2/Mw*(th.f_vvv-2*f_Tv/f_TT*th.f_Tvv+f_Tv^2/f_TT^2*th.f_TTv);
    J = [dhdT, dhdv, Mw*u;-th.f_TT, -th.f_Tv,0;dc2dT,dc2dv,-2*u];
end