function plotarticle(x, t, sol, PAR)
% Plotting the results
disp('Plotting Results')

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

R = 8.315; % Gas constant  J*mol^-1*K^-1;
F = 96500; % Faraday's constant (C/mol);
psifac = R*T/F; % Standard potential RT/F(V);

c_Ko = sol(:,:,1);
c_Ki = sol(:,:,2);
c_Nao = sol(:,:,3);
c_Nai = sol(:,:,4);
c_Clo = sol(:,:,5);
c_Cli = sol(:,:,6);
v_mo = sol(:,:,7);
v_mi = sol(:,:,8);


for j = 1:length(t)
%Derivative of c_Ko v/ at all t given by duoutdx1 (V/m);
[uout1(j,:),duoutdx1(j,:)] = pdeval(0,x,sol(j,:,1),x);
[uout2(j,:),duoutdx2(j,:)] = pdeval(0,x,sol(j,:,2),x);
[uout3(j,:),duoutdx3(j,:)] = pdeval(0,x,sol(j,:,3),x);
[uout4(j,:),duoutdx4(j,:)] = pdeval(0,x,sol(j,:,4),x);
[uout5(j,:),duoutdx5(j,:)] = pdeval(0,x,sol(j,:,5),x);
[uout6(j,:),duoutdx6(j,:)] = pdeval(0,x,sol(j,:,6),x);
[uout7(j,:),duoutdx7(j,:)] = pdeval(0,x,sol(j,:,7),x);
end

% Input/Output-vectors
Kinput = zeros(size(sol(:,:,1)));
Koutput = zeros(size(sol(:,:,1)));
I_tstart = min(find(t>=tstart));
I_tstop = max(find(t<tstop));
I_xstart = min(find(x>=xstart));
I_xstop = max(find(x<=xstop));
Kinput(I_tstart:I_tstop, I_xstart:I_xstop) = Imax*10^6; %Convert to micromol/(m^2s)
Koutput = Kdec*(c_Ko - c_Ko0)*10^6; %Convert to micromol/(m^2s)

%Flux densities due to diffusion (mol/(m^2*s));
j_KoD = -(D_K/lambda_o^2)*duoutdx1;
j_KiD = -(D_K/lambda_i^2)*duoutdx2;
j_NaoD = -(D_Na/lambda_o^2)*duoutdx3;
j_NaiD = -(D_Na/lambda_i^2)*duoutdx4;
j_CloD = -(D_Cl/lambda_o^2)*duoutdx5;
j_CliD = -(D_Cl/lambda_i^2)*duoutdx6;

%Current densities due to diffusion (A/m^2);
i_odiff = F*(j_KoD + j_NaoD - j_CloD);
i_idiff = F*(j_KiD + j_NaiD - j_CliD);

%Resistivities (ohm*m)
r_o = psifac*lambda_o^2./(F*(D_Na*c_Nao+D_K*c_Ko+D_Cl*c_Clo));
r_i = psifac*lambda_i^2./(F*(D_Na*c_Nai+D_K*c_Ki+D_Cl*c_Cli));

% Gradients of v_o og v_i (V/m)
dv_odx = (-duoutdx7 + r_i.*i_idiff+r_i*(A_o/A_i).*i_odiff).*r_o./(r_o+r_i*A_o/A_i);
dv_idx = (duoutdx7 + r_o.*(A_i/A_o).*i_idiff+r_o.*i_odiff).*(r_i./(r_i+r_o*A_i/A_o));

%Flux densities due to electrical migration(mol/(m^2*s));
j_KoV = -(D_K/lambda_o^2).*(1/psifac).*c_Ko.*dv_odx;
j_KiV = -(D_K/lambda_i^2).*(1/psifac).*c_Ki.*dv_idx;
j_NaoV = -(D_Na/lambda_o^2).*(1/psifac).*c_Nao.*dv_odx;
j_NaiV = -(D_Na/lambda_i^2).*(1/psifac).*c_Nai.*dv_idx;
j_CloV = (D_Cl/lambda_o^2).*(1/psifac).*c_Clo.*dv_odx;
j_CliV = (D_Cl/lambda_i^2).*(1/psifac).*c_Cli.*dv_idx;

%Current densities due to electrical migration (A/m^2);
i_ofield = F*(j_KoV + j_NaoV - j_CloV);
i_ifield = F*(j_KiV + j_NaiV - j_CliV);

% Total axial current densities (A/m^2);
i_i = i_idiff + i_ifield;
i_o = i_odiff + i_ofield;

%%% Membrane mechanisms
P = NaKpumprate(c_Ko, c_Nai); % Na/K pumprate (mol/(m^2*s));
f_Kir = Kirf(c_Ko,c_Ki,c_Ko0,c_Ki0, v_mo, v_m0, psifac); % Kir-factor (unitless)

% Membrane-flux densities (mol/(m^2*s));
j_KM = (g_K0.*f_Kir/F).*(v_mo-psifac.*log(c_Ko./c_Ki)) - 2*P;
j_NaM = (g_Na/F).*(v_mo-psifac.*log(c_Nao./c_Nai)) + 3*P;
j_ClM = -(g_Cl/F).*(v_mo+psifac.*log(c_Clo./c_Cli));

%Membrane current density (A/m^2);
i_m = F*(j_KM + j_NaM - j_ClM);

%% Nernst potentials (V)
eNa = psifac*log(c_Nao./c_Nai); 
eK = psifac*log(c_Ko./c_Ki);
eCl = -psifac*log(c_Clo./c_Cli);



% Indices for plotting
ItSS = max(find(t<tstop)); % time point while in steady state
Ixiz = 1; % x-value in input zone (x = 0)

%%% Some colors
ColK = [0.2 0.1 0.8]; % Blue
ColNa = [0.4 0.8 0.4]; %Green
ColCl = [0.8 0.2 0.7]; %Magenta
Cole = [1 0.7 0]; %Yellow
ColV = [0 0 0]; %Black
ColIN = [0 0 0]; %Black
Colr = [0 0 0]; % Black

x = x*1000; % convert to mm.


%% FIGURE 1: CONCENTRATION PROFILES
figure('position', [60 90 800 750], ...
    'name', 'Dynamics and Steady State profiles (Fig. 4 in Halnes et al. 2013)')

subplot(4,2,1) % Input
hold on;
plot(t, Kinput(:,Ixiz), 'Color', ColK, 'LineStyle', '-', 'LineWidth',2.0);
plot(t, -Koutput(:,Ixiz), 'Color', ColK, 'LineStyle', '--', 'LineWidth',2.0);
ylabel('j_K (\mu mol m^{-2}s^{-1})', 'FontSize', 12);
axis([0 600 -0.25 0.57]) 

subplot(4,2,3) % Extracellular concentrations in input zone (t)
c_Qo = c_Ko + c_Nao - c_Clo;
c_Qo0 = c_Ko0 + c_Nao0 - c_Clo0;
hold on
h=plot(t,c_Ko(:,Ixiz)-c_Ko0, 'Color', ColK);
set(h,'LineWidth',1.5);
h = plot(t, c_Nao(:,Ixiz)-c_Nao0, 'Color', ColNa);
set(h,'LineWidth',1.5);
h = plot(t, c_Clo(:,Ixiz)-c_Clo0, 'Color', ColCl);
set(h,'LineWidth',1.5);
h = plot(t, c_Qo(:,Ixiz) - c_Qo0, 'Color', Cole); 
set(h,'LineWidth',1.5);
ylabel('\Delta [k]_{E} (mM)', 'fontsize', 12);
axis([0 600 -35 10]);

subplot(4,2,5)%% Intracellular concentrations in input zone (t):
hold on
c_Qi = c_Ki + c_Nai - c_Cli;
c_Qi0 = c_Ki0 + c_Nai0 - c_Cli0;
hold on
h=plot(t,c_Ki(:,Ixiz)-c_Ki0, 'Color', ColK);
set(h,'LineWidth',1.5);
h = plot(t, c_Nai(:,Ixiz)-c_Nai0, 'Color', ColNa);
set(h,'LineWidth',1.5);
h = plot(t, c_Cli(:,Ixiz)-c_Cli0, 'Color', ColCl);
set(h,'LineWidth',1.5);
h = plot(t, c_Qi(:,Ixiz) - c_Qi0, 'Color', Cole); 
set(h,'LineWidth',1.5);
ylabel('\Delta [k]_{I} (mM)', 'fontsize', 12);
axis([0 600 -8 15])

subplot(4,2,7) % Voltage in input zone (t)
hold on;
h = plot(t, 1000*v_mo(:,[Ixiz])','Color', ColV);
set(h,'LineWidth',1.5);
xlabel('t (s)', 'FontSize', 12);
ylabel('v_M(mV)', 'FontSize', 12);
axis([0 600 -85 -55])

subplot(4,2,2) % Input/output during SS(x)
hold on
plot(x, Kinput(ItSS,:), 'Color', ColK, 'LineWidth', 1.5);
plot(x, -Koutput(ItSS,:),'Color', ColK, 'LineStyle', '--' , 'LineWidth', 1.5);
h = legend('j_K^{in}', 'j_K^{out}');
set(h, 'fontsize', 12)
axis([0 0.3 -0.25 0.57]);


subplot(4,2,4) % Extracellular concentrations during SS (x)
hold on
c_Qo = c_Ko + c_Nao - c_Clo;
c_Qo0 = c_Ko0 + c_Nao0 - c_Clo0;
h = plot(x,c_Ko(ItSS,:)-c_Ko0, 'Color', ColK);
set(h,'LineWidth',1.5);
h = plot(x, c_Nao(ItSS,:)-c_Nao0, 'Color', ColNa);
set(h,'LineWidth',1.5);
h = plot(x, c_Clo(ItSS,:)-c_Clo0, 'Color', ColCl);
set(h,'LineWidth',1.5);
h = plot(x, c_Qo(ItSS,:) - c_Qo0, 'Color', Cole);
set(h,'LineWidth',1.5);
axis([0 0.3 -35 10])

subplot(4,2,6) % Intracellular concentrations during SS (x)
c_Qi = c_Ki + c_Nai - c_Cli;
c_Qi0 = c_Ki0 + c_Nai0 - c_Cli0;
hold on
h = plot(x,c_Ki(ItSS,:)-c_Ki0, 'Color', ColK);
set(h,'LineWidth',1.5);
h = plot(x, c_Nai(ItSS,:)-c_Nai0, 'Color', ColNa);
set(h,'LineWidth',1.5);
h = plot(x, c_Cli(ItSS,:)-c_Cli0, 'Color', ColCl);
set(h,'LineWidth',1.5);
h = plot(x, c_Qi(ItSS,:) - c_Qi0, 'Color', Cole);
set(h,'LineWidth',1.5);
h = legend('K^+','Na^+', 'Cl^-', 'e^+');
set(h, 'fontsize', 12)
axis([0 0.3 -8 15])


subplot(4,2,8) % Voltage during SS (x)
hold on
h = plot(x, 1000*v_mo(ItSS,:)','Color', ColV);
set(h,'LineWidth',1.5);
xlabel('x (mm)', 'FontSize', 12);
axis([0 0.3 -85 -55])


%% FIGURE: MATTER CYCLING DURING SS
figure('position', [60 90 800 750], ...
    'name', 'Matter transports (flux densities) during steady state (Fig. 5 in Halnes et al. 2013)')

subplot(3,2,1) % Input/Output flux densities
hold on
set(gca, 'ColorOrder', [ColK;ColNa;ColCl;Cole]);
plot(x, Kinput(ItSS,:)-Koutput(ItSS,:), 'Color', ColK, 'LineWidth', 1.5);
plot(x, -Kinput(ItSS,:)+Koutput(ItSS,:), 'Color', ColNa, 'LineWidth', 1.5);
ylabel('j_k^{in}-j_k^{out}', 'fontsize', 12);
axis([0 0.3 -0.4 0.4]);

subplot(3,2,2) % Membrane flux densities
hold on
set(gca, 'ColorOrder', [ColK;ColNa;ColCl;Cole]);
j_QM = j_KM + j_NaM - j_ClM;
h = plot(x, 1e6*[j_KM(ItSS,:)', j_NaM(ItSS,:)', j_ClM(ItSS,:)', j_QM(ItSS,:)']);
set(h, 'LineWidth', 1.5);
ylabel('j_{kM}', 'fontsize', 12);
h = legend('K^+','Na^+', 'Cl^-', 'e^+');
set(h, 'fontsize', 10)
AX = axis; AX(2) = 0.3; axis(AX);
axis([0 0.3 -0.3 0.12]);


subplot(3,2,3)%% Extacellular field-flux densities:
hold on
j_QoV = (j_KoV + j_NaoV - j_CloV);
set(gca, 'ColorOrder', [ColK;ColNa;ColCl;Cole]);
h = plot(x, A_o*1e6*[j_KoV(ItSS,:)', j_NaoV(ItSS,:)', j_CloV(ItSS,:)', j_QoV(ItSS,:)']);
set(h, 'LineWidth', 1.5);
ylabel('a_Ej_{kE}^f', 'fontsize', 12);
axis([0 0.3 -100 55]);

subplot(3,2,4) % Intracellular field-flux densities
hold on
set(gca, 'ColorOrder', [ColK;ColNa;ColCl;Cole]);
j_QiV = (j_KiV + j_NaiV - j_CliV);
h = plot(x, A_i*1e6*[j_KiV(ItSS,:)', j_NaiV(ItSS,:)', j_CliV(ItSS,:)', j_QiV(ItSS,:)']);
set(h, 'LineWidth', 1.5);
ylabel('a_Ij_{kI}^f', 'fontsize', 12);
axis([0 0.3 -10 70]);

subplot(3,2,5) % Extracellular diffusion-flux densities
hold on
set(gca, 'ColorOrder', [ColK;ColNa;ColCl;Cole]);
j_QoD = (j_KoD + j_NaoD - j_CloD);
h = plot(x, A_o*1e6*[j_KoD(ItSS,:)', j_NaoD(ItSS,:)', j_CloD(ItSS,:)', j_QoD(ItSS,:)']);
set(h, 'LineWidth', 1.5);
ylabel('a_Ej_{kE}^d', 'fontsize', 12);
axis([0 0.3 -53 27]);
xlabel('x (mm)', 'fontsize', 12);

subplot(3,2,6) % Intracellular diffusion-flux densities
hold on
set(gca, 'ColorOrder', [ColK;ColNa;ColCl;Cole]);
j_QiD = j_KiD + j_NaiD - j_CliD;
h = plot(x, A_i*1e6*[j_KiD(ItSS,:)', j_NaiD(ItSS,:)', j_CliD(ItSS,:)', j_QiD(ItSS,:)']);
set(h, 'LineWidth', 1.5);
ylabel('a_Ij_{kI}^d', 'fontsize', 12);
axis([0 0.3 -3 9]);
xlabel('x (mm)', 'fontsize', 12);


%% FIGURE: Potassium uptake
figure('position', [60 90 400 750], ...
    'name', 'Potassium uptake/release during steady state (Fig. 6 in Halnes et al. 2013)')
AbsKirflux = (g_K0.*f_Kir/F).*(v_mo-psifac.*log(c_Ko./c_Ki))
AbsPumpflux =  2*P;

subplot(2, 1, 1)
hold on
h = plot(x, 1000*v_mo(ItSS,:)','k-');
set(h,'LineWidth',1.5);
h = plot(x, 1000*eK(ItSS,:)','b-');
set(h,'LineWidth',1.5);
xlabel('x (mm)', 'fontsize', 12);
ylabel('mV', 'FontSize', 12);
h = legend('v_M','e_K');
axis([0 0.3 -90 -55]);

subplot(2,1, 2)
hold on;
plot(x, 1e6*AbsKirflux(ItSS,:)', 'b--', 'LineWidth',2.0);
plot(x, 1e6*AbsPumpflux(ItSS,:)', 'b-', 'LineWidth',2.0);
xlabel('x (mm)', 'fontsize', 12);
ylabel('j (\mumol/(m^2 s))', 'fontsize', 12);
h = legend('|j_{Kir}|','|j_{Na^+/K^+ exchanger}|');
axis([0 0.3 0.55 1.1]);




