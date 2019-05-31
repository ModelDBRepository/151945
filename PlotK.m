function PlotK(x, t, sol, PAR)
% Plotting the results

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



%%%% PLOTS:

% Defining locations for plots
Itstart = min(find(t>tstart));
Itend = max(find(t<tstop));
Itmid = round((Itstart+Itend)/2);
Ixstart = min(find(x>xstart));
Ixend = max(find(x<xstop));
Ixmid = round((Ixstart+Ixend)/2);
Ixhalf = round(length(x)/2);
It = [Itstart-1, Itstart, Itmid, Itend, length(t)];
Ix = [1 Ixstart Ixmid, Ixend, length(x)];
Ixright = Ix(end);

%% Nernst potentials (V)
eNa = psifac*log(c_Nao./c_Nai); 
eK = psifac*log(c_Ko./c_Ki);
eCl = -psifac*log(c_Clo./c_Cli);


%% Input & Output
Ixmid = round((Ixstart+Ixend)/2);
Ixhalf = round(length(x)/2);
Ixright = round(length(x));

xA = round((Ixstart+Ixend)/2); % Middle of input zone
xB = round(length(x)/8); % Just outside input zone
xC = round(length(x)/2);  % Astrocyte midpoint

ton = max(find(t<102)); % Just after input onset
toff = max(find(t<402)); % just after input offset
tSS = max(find(t<300)); % some point during steady state
tAF = max(find(t<500)); % a while after stimulus was turned off


%%% Some colors
ColK = [0.2 0.1 0.8]; % Blue
ColNa = [0.4 0.8 0.4]; %Green
ColCl = [0.8 0.2 0.7]; %Magenta
Cole = [1 0.7 0]; %Yellow
ColV = [0 0 0]; %Black
ColIN = [0 0 0]; %Black
Colr = [0 0 0]; % Black


%% FIGURE 1: K
Kirflux = (g_K0.*f_Kir/F).*(v_mo-psifac.*log(c_Ko./c_Ki))
Pumpflux =  2*P;

figure
subplot(2, 1, 1)
hold on
% Voltage during SS (x)
h = plot(1000*x, 1000*v_mo(tSS,:)','k-');
set(h,'LineWidth',1.5);
h = plot(1000*x, 1000*eK(tSS,:)','b-');
set(h,'LineWidth',1.5);
xlabel('x (mm)', 'fontsize', 12);
ylabel('mV', 'FontSize', 12);
h = legend('v_M','e_K');
axis([0 0.3 -90 -55]);

subplot(2,1, 2)
hold on;
plot(1000*x, 1e6*Kirflux(tSS,:)', 'b--', 'LineWidth',2.0);
plot(1000*x, 1e6*Pumpflux(tSS,:)', 'b-', 'LineWidth',2.0);
%plot([0 0.3], [9.45e-1 9.45e-1], 'k--')
xlabel('x (mm)', 'fontsize', 12);
ylabel('j (\mumol/(m^2 s))', 'fontsize', 12);
h = legend('j_{Kir}','j_{Na^+/K^+ exchanger}');
axis([0 0.3 0.55 1.1]);



if 0
subplot(5,2,3) % Extracellular concentrations in input zone (t)
c_Qo = c_Ko + c_Nao - c_Clo;
c_Qo0 = c_Ko0 + c_Nao0 - c_Clo0;
hold on
h=plot(t,c_Ko(:,Ixmid)-c_Ko0, 'Color', ColK);
set(h,'LineWidth',1.5);
h = plot(t, c_Nao(:,Ixmid)-c_Nao0, 'Color', ColNa);
set(h,'LineWidth',1.5);
h = plot(t, c_Clo(:,Ixmid)-c_Clo0, 'Color', ColCl);
set(h,'LineWidth',1.5);
h = plot(t, c_Qo(:,Ixmid) - c_Qo0, 'Color', Cole); 
set(h,'LineWidth',1.5);
%xlabel('t (s)', 'FontSize', 12);
ylabel('\Delta [k]_{E} (mM)', 'fontsize', 12);
%Title('B', 'fontsize', 14);
%%%Title('c_{kE}(t,l/20) - c_{kE}_0', 'fontsize', 20);
%h = legend('K^+','Na^+', 'Cl^-', 'e^+');
%set(h, 'fontsize', 14)
axis([0 600 -35 10]);
%set(gca, 'FontSize', 12)


subplot(5,2,5)%% Intracellular concentrations in input zone (t):
hold on
c_Qi = c_Ki + c_Nai - c_Cli;
c_Qi0 = c_Ki0 + c_Nai0 - c_Cli0;
hold on
h=plot(t,c_Ki(:,Ixmid)-c_Ki0, 'Color', ColK);
set(h,'LineWidth',1.5);
h = plot(t, c_Nai(:,Ixmid)-c_Nai0, 'Color', ColNa);
set(h,'LineWidth',1.5);
h = plot(t, c_Cli(:,Ixmid)-c_Cli0, 'Color', ColCl);
set(h,'LineWidth',1.5);
h = plot(t, c_Qi(:,Ixmid) - c_Qi0, 'Color', Cole); 
set(h,'LineWidth',1.5);
%xlabel('t (s)', 'FontSize', 12);
%ylabel('mM', 'fontsize', 16);
ylabel('\Delta [k]_{I} (mM)', 'fontsize', 12);
%Title('C', 'fontsize', 14);
axis([0 600 -8 15])
%set(gca, 'FontSize', 16)


subplot(5,2,7) % Voltage in input zone (t)
hold on;
h = plot(t, 1000*v_mo(:,[xA])','Color', ColV);
%, t, 1000*eK(:,[xA])', 'b', ...
%t, 1000*eNa(:,[xA])', 'g', t, 1000*eCl(:,[xA])', 'r'); 
set(h,'LineWidth',1.5);
%xlabel('t (s)', 'FontSize', 12);
ylabel('v_M(mV)', 'FontSize', 12);
%h = legend('V_m(t,l/20)');
%set(h, 'fontsize', 16)
%Title('D', 'FontSize', 14);
axis([0 600 -85 -55])
%set(gca, 'FontSize', 16)


subplot(5,2,9) % resistivities in input zone (t)
hold on;
%h = plot(x, OutsQ(tSS,:)- (OutsQ(1,1)), 'Color', black, x, InsQ(tSS,:) - (InsQ(1,1)), 'Color', grey);
h = plot(t, r_o(:,Ixmid)'./r_o(1,1),'Color', Colr)
set(h, 'Linewidth', 1.5);
h = plot(t, r_i(:,Ixmid)'./r_i(1,1), 'Color', Colr);
set(h, 'Linewidth', 1.5, 'LineStyle','--');
xlabel('t (s)', 'FontSize', 12);
ylabel('r_n/r_{n0}', 'FontSize', 12);
%ylabel('a_n(\rho_n - \rho_{n0}) (C/m^3)', 'FontSize', 16);
%legend('n=E', 'n=I');
%Title('E', 'FontSize', 14)
%set(gca, 'FontSize', 16)
axis([0 600 0.85 1.25]);


subplot(5,2,2) % Input/output during SS(x)
hold on
plot(x, Kinput(tSS,:), 'Color', ColIN, 'LineWidth', 1.5);
plot(x, -Koutput(tSS,:),'Color', ColIN, 'LineStyle', '--' , 'LineWidth', 1.5);
%xlabel('x (mm)', 'FontSize', 12);
%ylabel('j^{in/out} ( \mu mol (m^{-2}s^{-1} )', 'FontSize', 12);
%Title('F', 'FontSize', 14);
h = legend('j^{in}', 'j^{out}');
set(h, 'fontsize', 12)
axis([0 0.3 -0.25 0.57]);
%set(gca, 'FontSize', 16)


subplot(5,2,4) % Extracellular concentrations during SS (x)
hold on
c_Qo = c_Ko + c_Nao - c_Clo;
c_Qo0 = c_Ko0 + c_Nao0 - c_Clo0;
h = plot(x,c_Ko(tSS,:)-c_Ko0, 'Color', ColK);
set(h,'LineWidth',1.5);
h = plot(x, c_Nao(tSS,:)-c_Nao0, 'Color', ColNa);
set(h,'LineWidth',1.5);
h = plot(x, c_Clo(tSS,:)-c_Clo0, 'Color', ColCl);
set(h,'LineWidth',1.5);
h = plot(x, c_Qo(tSS,:) - c_Qo0, 'Color', Cole);
set(h,'LineWidth',1.5);
%xlabel('x (mm)', 'FontSize', 12);
%ylabel('c_{kE} - c_{kE0} (mM)', 'fontsize', 12);
%Title('G', 'fontsize', 14);
%AX = axis; AX(2) = 0.3; axis(AX);
axis([0 0.3 -35 10])
%set(gca, 'FontSize', 16)

subplot(5,2,6) % Intracellular concentrations during SS (x)
c_Qi = c_Ki + c_Nai - c_Cli;
c_Qi0 = c_Ki0 + c_Nai0 - c_Cli0;
hold on
h = plot(x,c_Ki(tSS,:)-c_Ki0, 'Color', ColK);
set(h,'LineWidth',1.5);
h = plot(x, c_Nai(tSS,:)-c_Nai0, 'Color', ColNa);
set(h,'LineWidth',1.5);
h = plot(x, c_Cli(tSS,:)-c_Cli0, 'Color', ColCl);
set(h,'LineWidth',1.5);
h = plot(x, c_Qi(tSS,:) - c_Qi0, 'Color', Cole);
set(h,'LineWidth',1.5);
%xlabel('x (mm)', 'FontSize', 12);
%ylabel('c_{kI} - c_{kI0} (mM)', 'fontsize', 12);
%Title('H', 'fontsize', 14);
h = legend('K^+','Na^+', 'Cl^-', 'e^+');
set(h, 'fontsize', 12)
axis([0 0.3 -8 15])
%set(gca, 'FontSize', 16)



subplot(5,2,10) % Resistivities during SS (x)
hold on
h = plot(x, r_o(tSS,:)./r_o(1,1),'Color', Colr)
set(h, 'LineWidth', 1.5);
h = plot(x, r_i(tSS,:)./r_i(1,1), 'Color', Colr);
set(h, 'LineWidth', 1.5, 'LineStyle', '--');
xlabel('x (mm)', 'FontSize', 12);
%ylabel('r_n/r_{n0}', 'FontSize', 12);
%Title('J', 'FontSize', 14);
h = legend('ECS','ICS')
set(h, 'FontSize', 12);
%set(gca, 'FontSize', 16)
axis([0 0.3 0.85 1.25]);


%% FIGURE 2: MATTER CYCLING DURING SS
figure

subplot(3,2,1)%% Extacellular field-flux densities:
hold on
j_QoV = (j_KoV + j_NaoV - j_CloV);
set(gca, 'ColorOrder', [ColK;ColNa;ColCl;Cole]);
h = plot(x, A_o*1e6*[j_KoV(tSS,:)', j_NaoV(tSS,:)', j_CloV(tSS,:)', j_QoV(tSS,:)']);
set(h, 'LineWidth', 1.5);
%xlabel('x (mm)', 'fontsize', 12);
ylabel('a_Ej_{kE}^f', 'fontsize', 12);
%title('A', 'fontsize', 14);
%legend('K','Na', 'Cl');
axis([0 0.3 -100 55]);
%set(gca, 'FontSize', 12)

subplot(3,2,3) % Extracellular diffusion-flux densities
hold on
set(gca, 'ColorOrder', [ColK;ColNa;ColCl;Cole]);
j_QoD = (j_KoD + j_NaoD - j_CloD);
h = plot(x, A_o*1e6*[j_KoD(tSS,:)', j_NaoD(tSS,:)', j_CloD(tSS,:)', j_QoD(tSS,:)']);
set(h, 'LineWidth', 1.5);
%xlabel('x (mm)', 'fontsize', 12);
ylabel('a_Ej_{kE}^d', 'fontsize', 12);
%title('B', 'fontsize', 14);
%legend('K','Na', 'Cl');
axis([0 0.3 -53 27]);
%set(gca, 'FontSize', 16)

subplot(3,2,5) % Input/Output flux densities
hold on
set(gca, 'ColorOrder', [ColK;ColNa;ColCl;Cole]);
plot(x, Kinput(tSS,:)-Koutput(tSS,:), 'Color', ColK, 'LineWidth', 1.5);
plot(x, -Kinput(tSS,:)+Koutput(tSS,:), 'Color', ColNa, 'LineWidth', 1.5);
xlabel('x (mm)', 'fontsize', 12);
ylabel('j_k^{in}-j_k^{out}', 'fontsize', 12);
%title('C', 'fontsize', 14);
%legend('J_KM','J_NaM', 'J_ClM');
axis([0 0.3 -0.4 0.4]);
%set(gca, 'FontSize', 16)

subplot(3,2,2) % Intracellular field-flux densities
hold on
set(gca, 'ColorOrder', [ColK;ColNa;ColCl;Cole]);
j_QiV = (j_KiV + j_NaiV - j_CliV);
h = plot(x, A_i*1e6*[j_KiV(tSS,:)', j_NaiV(tSS,:)', j_CliV(tSS,:)', j_QiV(tSS,:)']);
set(h, 'LineWidth', 1.5);
%xlabel('x (mm)', 'fontsize', 12);
ylabel('a_Ij_{kI}^f', 'fontsize', 12);
%title('D', 'fontsize', 14);
axis([0 0.3 -10 70]);
%set(gca, 'FontSize', 12)

subplot(3,2,4) % Intracellular diffusion-flux densities
hold on
set(gca, 'ColorOrder', [ColK;ColNa;ColCl;Cole]);
j_QiD = j_KiD + j_NaiD - j_CliD;
h = plot(x, A_i*1e6*[j_KiD(tSS,:)', j_NaiD(tSS,:)', j_CliD(tSS,:)', j_QiD(tSS,:)']);
set(h, 'LineWidth', 1.5);
%xlabel('x (mm)', 'fontsize', 12);
ylabel('a_Ij_{kI}^d', 'fontsize', 12);
%title('E', 'fontsize', 14);
%legend('K','Na', 'Cl');
axis([0 0.3 -3 9]);
%set(gca, 'FontSize', 16)

end