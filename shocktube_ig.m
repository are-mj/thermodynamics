function [T2,p2,u2,T3,T5,p5,M1] ...
  = shocktube_ig(T1,p1,Mw1,gamma1,T4,p4,Mw4,gamma4)
% Shock tube model for a perfect gas
% Shock tube zones:
%    1: Undisturbed driven gas 
%    2: Shocked, driven gas
%    3: Driver gas between rarefaction fan and contact surface
%    4: Undisturbed driver gas
%    5: Reflected shock
% 
% References:
%  Kerwin, J.M.: Design of a Shock Tube for Jet Noise Research
%    MSc thesis, MIT, 1996. https://core.ac.uk/download/pdf/4414075.pdf
%  Landau, L.D. and Lifschitz, E.M.: Fluid Mechanics, 
%    Pergamon Press, 1959 (L&L)

  if nargin < 7
    Mw4 = Mw1;
    gamma4 = gamma1;
  end
  if nargin < 6
    T4 = T1;
  end
  R = 8314.3714;  % Universal gas constant
  a1 = sqrt(gamma1*R*T1/Mw1);  % Speed of sound
  a4 = sqrt(gamma4*R*T4/Mw4);
    
  ex = -2*gamma4/(gamma4-1);
  % Eq 2.1:
  fun = @(p2) p4-p2*(1-(gamma4-1)*a1/a4*(p2/p1-1)/ ...
    sqrt(2*gamma1*(2*gamma1+(gamma1+1)*(p2/p1-1))))^ex;
  % fun has a singularity,  make sure upperlimit is below this:
  r = roots([((gamma4-1)*a1/a4)^2,-2*gamma1*(gamma1+1),-4*gamma1^2]);
  upperlimit = p1*max(r);
  p2 = fzero(fun,[p1,upperlimit]);
  u2 = 2*a4/(gamma4-1)*(1-(p2/p4)^((gamma4-1)/2/gamma4)); % Eq. A8
  T3 = T4*(p2/p4)^((gamma4-1)/gamma4);
%   u3 = u2;
  M1 = sqrt((p2/p1*(gamma1+1)+(gamma1-1))/(2*gamma1));  % L&L (85.8) 
  T2 = T1*(2*gamma1*M1^2-(gamma1-1))*((gamma1-1)*M1^2 + 2) ...
    /(gamma1+1)^2/M1^2; % L&L (85.9)
  
  % Reflected shock (zone 5)
  % (L&L, Exercise 1, Paragraph 94)
  p5 = p2*((3*gamma1-1)*p2 - (gamma1-1)*p1) ...
    /((gamma1-1)*p2 + (gamma1+1)*p1);
  % L&L (85.2):
  T5 = T2*p5/p2*((gamma1+1)*p2 + (gamma1-1)*p5) ...
    /((gamma1-1)*p2 + (gamma1+1)*p5);
  % Shock speed ( from conservation of moles):
  % v5 = R*T5/p5;
  % W5 = u2*v5/(v2-v5);
end

