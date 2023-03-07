% Simulation of the filling of a gas tank from another tank at higher
% pressure
  th = thermo('H2');   % Create thermodynamic object
  % Example:
  Ta = 273.15+15;
  T_up = Ta;
  p_up = 350e5;
  T_down = Ta;
  p_down = 20e5;
  w = 0.2;  % kg/s
  [t,Z,pA,pB] = gas_tank_system(th,T_up,p_up,T_down,p_down,Ta,w);
  
  TA = Z(:,1);
  NA = Z(:,2);
  TB = Z(:,3);
  NB = Z(:,4);
  figure;
  subplot(311);
  plot(t,[TA,TB]-274.15)
  ax = axis; ax(3:4) = [min([TA;TB]),max([TA;TB])]-273.15;
  axis(ax);
  title Temperatures; ylabel \circC
  subplot(312);
  plot(t,[pA,pB]*1e-5);
  title Pressures;ylabel bar
  subplot(313);
  plot(t,[NA,NB]*th.Mw);
  title('Gas mass');ylabel kg
  ax = axis; ax(4) = ceil(max(NA)*th.Mw/25)*25;
  axis(ax);
  legend('Upstream tank','Downstream tank','location','west')
  xlabel seconds 