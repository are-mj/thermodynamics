% Demo of shocktube.m
% Comparing rigorous thermodynamics with perfect gas

R    = 8314.46261815324;     % Universal gas constant (J/(kmol K))
% Section1:   Low-pressure section
th1 = thermo('Air');   % Rigorous thermodynaics
Mw1 = 28.9586;
gamma1 = 1.4;         % Perfect gas
T1 = 300;
p1 = 1e5;

% Section 4: High-pressure section
th4 = thermo('H2');
Mw4 = 2.0159;
gamma4 = 1.4;
T4 = 300;
p4 = [1:9,10:10:700]'*1e5;

% allocate space
n = length(p4);
T2 = zeros(n,1);
T3 = zeros(n,1);
u2 = zeros(n,1);
p2 = zeros(n,1);
c2 = zeros(n,1);
c3 = zeros(n,1);
T5 = zeros(n,1);
p5 = zeros(n,1);
T2_ig = zeros(n,1);
T3_ig = zeros(n,1);
u2_ig = zeros(n,1);
p2_ig = zeros(n,1);
c2_ig = zeros(n,1);
c3_ig = zeros(n,1);
T5_ig = zeros(n,1);
p5_ig = zeros(n,1);

% Solve for range of driver pressures
tic;
for i = 1:n
  [T2(i),p2(i),u2(i),T3(i),T5(i),p5(i)] = shocktube(T1,p1,th1,T4,p4(i),th4);
  th1.Tpcalc(T2(i),p2(i));
  c2(i) = th1.c;
  th4.Tpcalc(T3(i),p2(i));
  c3(i) = th4.c;
  [T2_ig(i),p2_ig(i),u2_ig(i),T3_ig(i),T5_ig(i),p5_ig(i)] ...
    = shocktube_ig(T1,p1,Mw1,gamma1,T4,p4(i),Mw4,gamma4);
  c2_ig(i) = sqrt(gamma1*R*T2(i)/Mw1);
  c3_ig(i) = sqrt(gamma4*R*T3(i)/Mw4);
end
toc

figure;
h = plot(p4*1e-5,[p2,p5,p2_ig,p5_ig]*1e-5);
for i = 1:2
  col = get(h(i),'color');
  set(h(i+2),'color',col,'linestyle','--');
end
grid on;
xlabel('Driver gas pressure (bar)')
ylabel('Pressures (bar)')
set(gca,'Gridalpha',0.3)
title(sprintf('Shock tube pressures. Driver gas: %s, Driven gas: %s',th4.species,th1.species))
legend('Shock','Reflected shock','Perfect gas','location','east')

figure;
h = semilogy(p4*1e-5,[T2,T3,T5,T2_ig,T3_ig,T5_ig]);
for i = 1:3
  col = get(h(i),'color');
  set(h(i+3),'color',col,'linestyle','--');
end
grid on;
xlabel('Driver gas pressure (bar)') 
ylabel('Temperature (K)')
set(gca,'ytick',[200,300,500,1000,2000,3000,5000])
set(gca,'ylim',[100,8000])
set(gca,'Gridalpha',0.3)
title(sprintf('Shock tube temperatures. Driver gas: %s, Driven gas: %s',th4.species,th1.species))
legend('Behind shock','Expanded driver gas','Behind reflected shock','Perfect gas','location','east')

%{
figure;
h = plot(p4*1e-5,[u2,c2,c3,u2_ig,c2_ig,c3_ig]);
for i = 1:3
  col = get(h(i),'color');
  set(h(i+3),'color',col,'linestyle','--');
end
title(sprintf('Velocities. Driver gas: %s, Driven gas: %s',th4.species,th1.species))
legend('Flow velocity','Zone 2 sped of sound','zone 3 speed of sound','Perfect gas','location','southeast')
%}