% Simulate the filling of a hydrogen tank with constant rate

  th = thermo('H2');      % Create thermodynamic object
  % Tank parameters:
  tank.wallmass = 500;    % kg
  tank.cp_wall = 100;     % J/kg/K 
  tank.V = 5;             % m3
  tank.Area = 2;          % m2
  tank.heatcoeff  = 0;    % W/m2/K  (Adiabatic tank)
  
  Ta = 273.15+15;   % ambient temoerature
  T0 = Ta;
  p0 = 20e5;
  T_up = Ta;
  p_up = 350e5;
  w = 0.2;  % Filling rate (kg/s)
  p_final = 300e5;
  [t,T,N,p] = fill_gas_tank(th,tank,T0,p0,T_up,p_up,Ta,w,p_final);

  figure;
  subplot(311);
  plot(t,T-274.15)
  ylabel \circC
  title('Filling of H_2 tank with constant rate')
  subplot(312);
  plot(t,p*1e-5);
  ylabel bar
  subplot(313);
  plot(t,N*th.Mw);
  ylabel kg  
  xlabel seconds  
  
  fprintf('Final temperature %6.1f\n',T(end)-273.15)
  fprintf('Temperature increase %6.1f\n',T(end)-T(1))