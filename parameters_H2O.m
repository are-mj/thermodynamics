function par = parameters_H2O
% Parameters for H2O to be used with helmholtz.m
% and with the 'thermo' thermodynamic object.  
% Sources: 
%    General thermodynamic properties: NIST webbook
%    Thermodynamic model parameters:  
%      W.Wagner, A.Pruss: J. Phys. Chem. Ref. Data 31, 387 (2002);

  par.species = 'H2O';
  par.casno = '7732185';
  % Data from NIST:
%   par.R    = 8314.46261815324;     % Universal gas constant (J/(kmol K))
%   par.Mw   = 18.0153;       % Molar mass (kg/kmol)
%   par.Tc   = 647;           % Critical temperature (K)
%   par.rhoc = 17.91;         % Critical density (NIST - kmol/m3)
%   par.vc   = 1/par.rhoc;    % Crtical molar volume (m3/kmol)
%   par.pc   = 220.64e5;      % Critical pressure (Pa)
%   par.Tt   = 273.16;        % Triple point temperature (K)
%   par.pt   = 611.655;       % Triple point pressure (Pa)
%   par.vlt  = 0.018019;      % Liquid molar volume at triple point (m3/kmol)
%   par.vvt  = 3710.98;       % Vapur phase molar volume at triple point
  % data from Wagner & Pruss (2002)
  par.R    = 8314.3714;     % Universl gas constant (J/(kmol K)
  par.Mw   = 18.015268;     % Molar mass (kg/kmol)
  par.Tc   = 647.096;       % Critical temperature (K)
  par.rhoc = 322/par.Mw;    % Critical density (kmol/m3)
  par.vc   = 1/par.rhoc;    % Crtical molar volume (m3/kmol)
  par.pc   = 220.64e5;      % Critical pressure (Pa)
  par.Tt   = 273.16;        % Triple point temperature (K)
  par.pt   = 611.655;       % Triple point pressure (Pa)
%   par.vlt  = 0.018019;      % Liquid molar volume at triple point (m3/kmol)
%   par.vvt  = 3710.98;       % Vapur phase molar volume at triple point
  
  % Ideal gas (Table 6.1):
  par.ig_a = [
    -8.32044648201
    6.6832105268    
    3.00632         
    0.012436  
    0.97315
    1.27950
    0.96956
    0.24873
    ];
  par.ig_a(1:2) = par.ig_a(1:2) + [4754.95;-2485690/par.Tc]/par.R;
  par.ig_b = [
    1.28728967
    3.53734222
    7.74073708
    9.24437796
    27.5075105];
  
  % Residual contribution  (Table 31)
  par.sections = [7,51,54,56];
  par.n = [
   0.12533547935523e-1
   0.78957634722828e1
  -0.87803203303561e1
   0.31802509345418
  -0.26145533859358
  -0.78199751687981e-2
   0.88089493102134e-2
  -0.66856572307965
   0.20433810950965
  -0.66212605039687e-4
  -0.19232721156002
  -0.25709043003438
   0.16074868486251
  -0.40092828925807e-1
   0.39343422603254e-6
  -0.75941377088144e-5
   0.56250979351888e-3
  -0.15608652257135e-4
   0.11537996422951e-8
   0.36582165144204e-6
  -0.13251180074668e-11
  -0.62639586912454e-9
  -0.10793600908932
   0.17611491008752e-1
   0.22132295167546
  -0.40247669763528
   0.58083399985759
   0.49969146990806e-2
  -0.31358700712549e-1
  -0.74315929710341
   0.47807329915480
   0.20527940895948e-1
  -0.13636435110343
   0.14180634400617e-1
   0.83326504880713e-2
  -0.29052336009585e-1
   0.38615085574206e-1
  -0.20393486513704e-1
  -0.16554050063734e-2
   0.19955571979541e-2
   0.15870308324157e-3
  -0.16388568342530e-4
   0.43613615723811e-1
   0.34994005463765e-1
  -0.76788197844621e-1
   0.22446277332006e-1
  -0.62689710414685e-4
  -0.55711118565645e-9
  -0.19905718354408
   0.31777497330738
  -0.11841182425981
  -0.31306260323435e2  
   0.31546140237781e2  
  -0.25213154341695e4  
  -0.14874640856724    
   0.31806110878444];  


  par.d = [
  1 
  1 
  1 
  2 
  2 
  3 
  4 
  1 
  1 
  1 
  2 
  2 
  3 
  4 
  4 
  5 
  7 
  9 
  10
  11
  13
  15
  1 
  2 
  2 
  2 
  3 
  4 
  4 
  4 
  5 
  6 
  6 
  7 
  9 
  9 
  9 
  9 
  9 
  10
  10
  12
  3 
  4 
  4 
  5 
  14
  3 
  6 
  6 
  6 
  3 
  3 
  3 
  ];

  par.a = [
  3.500 
  3.500 ];

  par.t = [
  -0.5
  0.875
  1   
  0.5 
  0.75
  0.375
  1   
  4   
  6   
  12  
  1   
  5   
  4   
  2   
  13  
  9   
  3   
  4   
  11  
  4   
  13  
  1   
  7   
  1   
  9   
  10  
  10  
  3   
  7   
  10  
  10  
  6   
  10  
  10  
  1   
  2   
  3   
  4   
  8   
  6   
  9   
  8   
  16  
  22  
  23  
  23  
  10  
  50  
  44  
  46  
  50 
  0
  1
  4
   ];

  par.b = [
  0.85 
  0.95 
  ];

  par.c = [
  1  
  1  
  1  
  1  
  1  
  1  
  1  
  1  
  1  
  1  
  1  
  1  
  1  
  1  
  1  
  2  
  2  
  2  
  2  
  2  
  2  
  2  
  2  
  2  
  2  
  2  
  2  
  2  
  2  
  2  
  2  
  2  
  2  
  2  
  2  
  3  
  3  
  3  
  3  
  4  
  6  
  6  
  6  
  6  
  ];    

  par.alpha = [20;20;20];
  par.beta = [150;150;250;0.3;0.3];
  par.gamma = [1.21;1.21;1.25];
  par.epsilon = ones(3,1);
  par.A = [0.32;0.32];
  par.B = [0.2;0.2];
  par.C = [28;32];
  par.D = [700;800];
  
  % Saturation pressure data (Table 8)   (ps = pc*exp(Tc/T*as*theta.^ase); 
  par.as = [-7.85951783,1.84408259,-11.7866497,22.6807411,-15.9618719,...
     1.80122502];
  par.ase = [1;1.5;3;3.5;4;7.5];
  % Saturated liquid volume  (vl = vc/(bs*theta.^bse)
  par.bs = [1,1.99274064,1.09965342,-0.510839303,-1.75493479,...
     -45.5170352,-6.74694450e5];
  par.bse = [0;1;2;5;16;43;110]/3;     
  par.cs = [2.03150240,2.68302940,5.38626492,17.2991605,44.7586581,...
      63.9201063];
  par.cse = [2;4;8;18;37;71]/6;    
  % Saturated liquid volume (vv = vc/(cs*theta.^cse)
  
  % Vapour/solid sublimation curve
  % Wagner, W. et al., J. Phys. Chem. Ref. Data 40.4 (2011)
  par.sublima = [-0.212144006e2,0.273203819e2,-0.610598130e1];
  par.sublime = [0.333333333e-2,0.120666667e1,0.170333333e1];
  
  % Liquid/solid melting curve  
  % Wagner, W. et al., J. Phys. Chem. Ref. Data 40.4 (2011)
  par.melta = [0.119539337e7,0.808183159e5,0.333826860e4];
  par.melte = [0.300000e1,0.257500e2,0.103750e3]; % Exponents  
end
