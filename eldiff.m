function [x, t, S, sol, PAR] = eldiff
%%%%% Astrocyte model
%%%%% Implemented by Geir Halnes, March 2013 

%%% PHYSICAL CONSTANTS %%%%%%%%%%%%%%%%%% Units: mol, m, V, s, C, ohm, A
R = 8.315; % Gas constant  J*mol^-1*K^-1;
F = 96500; % Faraday's constant (C/mol);
T = 298; % Temperature, K;
psifac = R*T/F; %Standard potential (V);

%%%% GEOMETRICAL PARAMETERS
% Notation here: o means outside (ECS), i means inside (intracellular)
l = 3.00e-4; %Length of cell (m);
A_i = 0.4; % Chen&Nicholson (Tissue volume fraction being astrocytes)
A_o = 0.2; %Chen&Nicholson (Tissue volume fraction being ECS)
O_m = A_i/5.00e-8; % (Astrocytic membrane area per tissue volume)

%%%% INPUT PARAMETERS
tstart = 100; tstop = 400; % Input starts and ends (s)
%tstart = 10; tstop = 45; % Used with compare_models.m
xstart = 0; 
xstop = l*0.1; % Region for input (m)

Imax = 5.5e-7; % Input amplitude K/Na-(exchange) (mol/(s m^2));
Kdec = 2.9e-8; % Output factor (m/s)

%%% DISCRETISATION
simt = 600; % Simulate for 600 s
xres = 1*100; x = linspace(0,l,xres); % 100 x-points
tres = 1*100;  t = linspace(0,simt,tres); % 100 t-points


%%% MEMBRANE PARAMETERS
g_Na = 1; % Na+ conductance (S*m^-2);
g_Cl=0.5; % Cl- concuctance (S*m^-2);
g_K0=16.96; % K+ conductance (S*m^-2);
C_m = 1.00e-2; % Membrane capacitance (Farad/m^2);

%%%% ELECTRODIFFUSION PARAMETERS (from Grodzinski-book)
D_K = 1.96e-9; %Diffusion coefficients (m^2/s);
D_Cl=2.03e-9;
D_Na=1.33e-9;
lambda_o = 1.6; % ECS tortuousity (1), C&N ;
lambda_i=3.2; % Intracellular tortuousity (1), C&N;

%%% STATIC CHARGES GIVING THE CORRECT CONCENTRATION/V_M-relation
u0 = CN34ic; % Loading the initial conditions
c_Ko0 = u0(1); c_Ki0 = u0(2); c_Nao0 = u0(3); c_Nai0 = u0(4);
c_Clo0 = u0(5); c_Cli0 = u0(6); v_m0 = u0(7);
X_iZ_i = -(F*A_i/O_m)*(c_Ki0 + c_Nai0 - c_Cli0) + C_m*v_m0; %Intracellular static charge (C/m^2)
X_oZ_o = -(F*A_o/O_m)*(c_Ko0 + c_Nao0 - c_Clo0) - C_m*v_m0; %ECS static charge (C/m^2)

%%%% Gather all parameters in vector named PAR
PAR = [l; A_i; A_o; O_m; ...
g_K0; g_Na; g_Cl; C_m; ...
c_Ko0; c_Ki0; c_Nao0; c_Nai0; c_Clo0; c_Cli0; ...
v_m0; X_iZ_i; X_oZ_o; ...
R; F; T; psifac; ...
D_K; D_Na; D_Cl; lambda_i; lambda_o;...
tstart; tstop; xstart; xstop; ...
Imax; Kdec;];

disp('starting simulation')
disp(['astrocyte length = ', num2str(l*1000), 'mm'])
disp(['input zone from x = ', num2str(xstart*1000), 'to', num2str(xstop*1000), ' mm'])
disp(['input applied from t = ', num2str(tstart), 's to t = ', num2str(tstop), ' s.'])
disp(['will simulate', num2str(simt), 's of activity'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Call pdepe-solver;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = odeset('RelTol',1e-4,'AbsTol',1e-8,'NormControl','on');
options = odeset('RelTol',1e-4,'AbsTol',1e-8,'NormControl','on', 'MaxStep', 0.1);

m = 0; %pdepe x-parameter;
sol = pdepe(m,@(x,t,u,DuDx)CN34pde(x,t,u,DuDx,PAR),...
    @(x)CN34ic,@(xl,ul,xr,ur,t)CN34bc(xl,ul,xr,ur,t,PAR),x,t,options);

% Store variables in the structure-variable S
S.c_Ko = sol(:,:,1);
S.c_Ki = sol(:,:,2);
S.c_Nao = sol(:,:,3);
S.c_Nai = sol(:,:,4);
S.c_Clo = sol(:,:,5);
S.c_Cli = sol(:,:,6);
S.v_mo = sol(:,:,7);
S.v_mi = sol(:,:,8);
disp('done!!');


function [c,f,s] = CN34pde(x,t,u,DuDx,PAR)
%CN34pde solves the pde-system: c Du/Dt = D/Dx[f(u,Du/Dx)] + s;

% Unwrapping parameters
l = PAR(1); A_i = PAR(2); A_o = PAR(3); O_m = PAR(4); 
g_K0 = PAR(5); g_Na = PAR(6); g_Cl = PAR(7); C_m = PAR(8);
c_Ko0 = PAR(9); c_Ki0 = PAR(10); c_Nao0 = PAR(11); c_Nai0 = PAR(12); 
c_Clo0 = PAR(13); c_Cli0 = PAR(14); 
v_m0 = PAR(15); X_iZ_i = PAR(16); X_oZ_o = PAR(17);
R = PAR(18); F = PAR(19); T = PAR(20); psifac = PAR(21);
D_K = PAR(22); D_Na = PAR(23); D_Cl = PAR(24); 
lambda_i = PAR(25); lambda_o = PAR(26);
tstart = PAR(27); tstop = PAR(28); xstart = PAR(29); xstop = PAR(30);
Imax = PAR(31); Kdec = PAR(32);

%%% Variable values (V, mol/m^3,A/m^2);
%%% o means outside (ECS), i means inside (intracellular)
u(find(u(1:6)<0)) = 0;
c_Ko = u(1); c_Ki = u(2); 
c_Nao = u(3); c_Nai = u(4); 
c_Clo = u(5); c_Cli = u(6);


%%% Algebraic method (define V_m in terms of ion concentrations)
v_mA = -(F*(A_o/O_m)/C_m*(c_Ko+c_Nao-c_Clo)+ X_oZ_o/C_m); %v_m derived from extracellular charge
v_mB = +(F*(A_i/O_m)/C_m*(c_Ki+c_Nai-c_Cli)+ X_iZ_i/C_m); %v_m defived from intracellular charge
v_m = v_mA; % Random choice (as v_mA = v_mB)
vmgrad_A = -F*(A_o/O_m)/C_m*(DuDx(1)+DuDx(3)-DuDx(5)); % v_m gradient
vmgrad_B = +F*(A_i/O_m)/C_m*(DuDx(2)+DuDx(4)-DuDx(6)); % v_m gradients
vmgrad = vmgrad_A; % Random choice (as v_mA = v_mB)

%%% Differential method (compute v_M from differential equation)
%    v_mA = u(7);
%    v_mB = u(8);
%    vmgrad_A = DuDx(7);
%    vmgrad_B = DuDx(8);
%    v_m = v_mA; % Random choice (Should have: v_mi = v_mo)
%    vmgrad = vmgrad_A;

%%% MEMBRANE MECHANISMS
% Kir-faktor for K+ conductance (unitless);
f_Kir = Kirf(c_Ko,c_Ki,c_Ko0,c_Ki0,v_m,v_m0,psifac);

% Na/K-pumperate;
P = NaKpumprate(c_Ko, c_Nai);

%%% INPUT:
if t>tstart & t<tstop & x>xstart & x<xstop
    ft = Imax - Kdec*(c_Ko-c_Ko0); % In the input zone, in the input time-window
else
    ft = - Kdec*(c_Ko-c_Ko0); % Decay towards resting concentrations
end

% Membrane flux densities (mol/(m^2*s));
j_KM = (g_K0*f_Kir/F)*(v_m-psifac*log(c_Ko/c_Ki)) - 2*P;
j_NaM = (g_Na/F)*(v_m-psifac*log(c_Nao/c_Nai)) + 3*P;
j_ClM = -(g_Cl/F)*(v_m+psifac*log(c_Clo/c_Cli));

% Membrane current density (A/m^2);
i_m = F*(j_KM+j_NaM-j_ClM);

% Flux densities due to diffusion (mol/(m^2*s));
j_KoD = -(D_K/lambda_o^2)*DuDx(1);
j_KiD = -(D_K/lambda_i^2)*DuDx(2);
j_NaoD = -(D_Na/lambda_o^2)*DuDx(3);
j_NaiD = -(D_Na/lambda_i^2)*DuDx(4);
j_CloD = -(D_Cl/lambda_o^2)*DuDx(5);
j_CliD = -(D_Cl/lambda_i^2)*DuDx(6);

% Current densities due to diffusion (A/m^2);
i_odiff = F*(j_KoD + j_NaoD - j_CloD);
i_idiff = F*(j_KiD + j_NaiD - j_CliD);

% Resistivities (Ohm*m)
r_o = psifac*lambda_o^2/(F*(D_Na*c_Nao+D_K*c_Ko+D_Cl*c_Clo));
r_i = psifac*lambda_i^2/(F*(D_Na*c_Nai+D_K*c_Ki+D_Cl*c_Cli));

%Gradients of v_o (ECS) and v_i (ICS), units (V/m);
dv_odx = (-vmgrad+r_i*i_idiff+r_i*(A_o/A_i)*i_odiff)*r_o/(r_o+r_i*A_o/A_i);
dv_idx = (vmgrad+r_o*(A_i/A_o)*i_idiff+r_o*i_odiff)*r_i/(r_i+r_o*A_i/A_o);

% Flux densities due to elecrical migration(mol/(m^2*s));
j_KoV = -(D_K/lambda_o^2)*(1/psifac)*c_Ko*dv_odx;
j_KiV = -(D_K/lambda_i^2)*(1/psifac)*c_Ki*dv_idx;
j_NaoV = -(D_Na/lambda_o^2)*(1/psifac)*c_Nao*dv_odx;
j_NaiV = -(D_Na/lambda_i^2)*(1/psifac)*c_Nai*dv_idx;
j_CloV = (D_Cl/lambda_o^2)*(1/psifac)*c_Clo*dv_odx;
j_CliV = (D_Cl/lambda_i^2)*(1/psifac)*c_Cli*dv_idx;

% Current densities due to electrical migration (A/m^2);
i_ofield = F*(j_KoV + j_NaoV - j_CloV);
i_ifield = F*(j_KiV + j_NaiV - j_CliV);

% Sanity checks:
if abs(i_ofield/(1/r_o*dv_odx))>1.01 | abs(i_ofield/(1/r_o*dv_odx))<0.99
disp('warning: field currents not consistent. Relative difference:')
i_ofield/(-1/r_o*dv_odx)
end
if abs(i_ifield/(1/r_i*dv_idx))>1.01 | abs(i_ifield/(1/r_i*dv_idx))<0.99
disp('warning: field currents not consistent. Relative difference:')
i_ifield/(-1/r_i*dv_idx)
end

% CN34pde solves the pde-system: c Du/Dt = D/Dx[f(u,Du/Dx)] + s ;
% Here we define c, f and s:
c = [1; 1; 1; 1; 1; 1; -(O_m/A_o)*C_m; +(O_m/A_i)*C_m];
s1 = (O_m/A_o)*(ft + j_KM);
s2 = -(O_m/A_i)*j_KM;
s3 = (O_m/A_o)*(-ft +j_NaM);
s4 = -(O_m/A_i)*j_NaM;
s5 = (O_m/A_o)*j_ClM;
s6 = -(O_m/A_i)*j_ClM;
s7 = (O_m/A_o)*i_m;
s8 = -(O_m/A_i)*i_m;
s=[s1; s2; s3; s4; s5; s6; s7; s8];
f1 = -j_KoD-j_KoV;
f2 = -j_KiD-j_KiV;
f3 = -j_NaoD-j_NaoV;
f4 = -j_NaiD-j_NaiV;
f5 = -j_CloD-j_CloV;
f6 = -j_CliD-j_CliV;
f7 = -i_odiff-i_ofield;
f8 = -i_idiff-i_ifield;
f = [f1; f2; f3; f4; f5; f6; f7; f8];

% Spit out some values once in a while during the simulation.
if rand>0.9999
disp(['simulated: ', num2str(t), ' s.'])
end


function u0 = CN34ic
% Initial conditions (mol/m^3 = mM);
% c_kn0 = Value from literature + correction (giving steady state)
c_Ko0=3.0+0.082;
c_Ki0=100.0 - 0.041;
c_Nao0=145.0 - 0.338;
c_Nai0=15.0 + 0.189;
c_Clo0=134 - 0.29;
c_Cli0=5 + 0.145;
v_m0 = -83.6e-3; % Membrane pot. (V);
u0 = [c_Ko0 c_Ki0 c_Nao0 c_Nai0 c_Clo0 c_Cli0, v_m0, v_m0]';


function [pl,ql,pr,qr] = CN34bc(xl,ul,xr,ur,t,param)
% Boundary conditions on the form: p + q f = 0;
% We have: f = [jKi, jKo, jNai, jNao, jCli, jClo]'
% Assuming zero flux/current through ends of cell; f = 0.
pl = [0; 0; 0; 0; 0; 0; 0; 0]; %left
pr = [0; 0; 0; 0; 0; 0; 0; 0]; % right
ql = [1; 1; 1; 1; 1; 1; 1; 1]; %left
qr = [1; 1; 1; 1; 1; 1; 1; 1]; % right

