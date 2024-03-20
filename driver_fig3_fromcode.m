% Run simulation of clinical PKPD model
clear all;

%% Options
do_biomarkers = 0;

%% set initial condition
CARTe_PB0 = 0; % CARTe in blood
CARTm_PB0 = 0; % CARTm in blood
CARTe_T0  = 0; % CARTe in tissue
CARTm_T0  = 0; % CARTm in tissue
Cplx0     = 0; % CAR-Target Complexes
Tumor0    = 1e8; %1e9; %2e9; %2.5e9; % tumor size

tmp = 12.1*2.5e9;
M0 = tmp/0.117;
B0 = 0.5;

IC = [CARTe_PB0;CARTm_PB0;CARTe_T0;CARTm_T0;Cplx0;Tumor0;M0;B0];

%% set parameters
p = set_params('PKPDclin_with_biomarkers');
p.Tumor0 = Tumor0;
p.M0 = M0;
p.B0 = B0;
% change params here

%p.Kkill_max = 0.05*p.Kkill_max; %0.1*p.Kkill_max; %0.06*p.Kkill_max;

[params, ~] = pars2vector(p, 0);


%% time span
t0 = 0;
tf = 360; 
tspan = [t0,tf];

%% CART dose
doseCART_tot = 50e6; %150e6; %450e6; %800e6; %450e6; %150e6;%50e6; %; % total number of cells in dose

%% ODE settings
options = odeset('RelTol',1.0e-12,'AbsTol',1e-12); % ode solver settings

%% Simulation
fprintf('treatment simulation \n')
IC(1) = doseCART_tot; % add dose into CARTe_PB
[t_treat, y_treat] = ode15s(@(t,y) modeqns_PKPD_withBM(t,y,params),...
                                    tspan, IC, options);
fprintf('sim finished \n')

fprintf('vehicle simulation \n')
IC(1) = 0;
[t_veh,y_veh] = ode15s(@(t,y) modeqns_PKPD_withBM(t,y,params),...
                                    tspan, IC, options);
fprintf('sim finished \n')


%% Make figures
fprintf('making figs \n')
% figure specs
lw = 4;
f.xlab = 16; f.ylab = 16; f.title = 18;
f.leg = 16; f.gca = 18;
cmap = parula(5);
c1 = cmap(1,:);
c2 = cmap(3,:);
ls1 = '-'; ls2 = '-';
cgraymap = gray(5);
cgray = cgraymap(3,:);
lwgray = 2; lsgray = '--';
labs = {'vehicle', 'treatment'};

%ymin = 0;
ymax = 1;

%% treatment only
figure(1)
nr = 2;
nc = 2;
clf;
subplot(nr,nc,1)
plot(t_treat,max(1e-16,y_treat(:,1)),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTe_{PB}')
%ylim([10^0, 10^8])
set(gca,'fontsize',f.gca,'YScale','log')
grid on

subplot(nr,nc,2)
plot(t_treat,y_treat(:,2),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTm_{PB}')
%ylim([10^0, 10^8])
set(gca,'fontsize',f.gca,'YScale','log')
grid on

subplot(nr,nc,3)
plot(t_treat,y_treat(:,3),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTe_{T}')
%ylim([10^0, 10^8])
set(gca,'fontsize',f.gca,'YScale','log')
grid on

subplot(nr,nc,4)
plot(t_treat,y_treat(:,4),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTm_{T}')
%ylim([10^0, 10^8])
set(gca,'fontsize',f.gca,'YScale','log')
grid on


temp = sprintf('Treatment dose %0.1e cells', doseCART_tot);
sgtitle(temp)


%% tumor
figure(2)
clf;
hold on
plot(t_veh,y_veh(:,6)/Tumor0,'linewidth',lw,'color',c1)
plot(t_treat,y_treat(:,6)/Tumor0,'linewidth',lw,'color',c2,'linestyle',':')
xlabel('t')
ylabel('Tumor/Tumor_0')
legend('Vehicle', 'Treatment')
set(gca,'fontsize',f.gca)
grid on

%% Blood transgene
TransC=0.002; % from model code
CARTe_PB = y_treat(:,1);
CARTm_PB = y_treat(:,2);
FinalCARTPB = TransC * (CARTe_PB + CARTm_PB);
figure(3)
clf
plot(t_treat, FinalCARTPB, 'linewidth',lw)
%plot(t_treat, FinalCARTPB, 'linewidth',lw,'color',c2)
xlabel('t (days)')
ylabel('Transgene copies/\mu g genomic DNA')
xlim([0,350])
ylim([10^-5, 10^7])
set(gca,'fontsize',f.gca, 'Yscale', 'log')
grid on

%% Cplx
figure(4)
clf;
plot(t_treat, y_treat(:,5),'linewidth',lw,'color',c2)
xlabel('t (days)')
ylabel('Cplx')
set(gca,'fontsize',f.gca)
grid on


%% Biomarkers
figure(5)
clf;
nr = 1; nc = 3;
subplot(nr,nc,1)
Mchange = (y_treat(:,7) - M0)/M0 * 100;
plot(t_treat, Mchange,'linewidth',lw,'color',c2)
xlabel('t (days)')
ylabel('M change')
xlim([0,60])
ylim([-100,100])
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,2)
Bchange = (y_treat(:,8) - B0)/B0 * 100;
plot(t_treat, Bchange,'linewidth',lw,'color',c2)
xlabel('t (days)')
ylabel('B change')
xlim([0,60])
ylim([-100,100])
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,3)
plot(t_treat,y_treat(:,6)/Tumor0,'linewidth',lw,'color',c2,'linestyle','-')
xlabel('t')
ylabel('Tumor/Tumor_0')
xlim([0,60])
set(gca,'fontsize',f.gca)
grid on

y=y_treat;
t=t_treat;
run('compute_alg_eqns.m')

%% save simulations options
save_sim = input('Do you want to save the simulation? (0 - no/1 - yes) ');
if save_sim
    notes = input('notes: ');
    fname = strcat('./simPKPD/', date, '_driver_fig3_fromcode', ...
                    '_dose-', num2str(doseCART_tot),...
                    '_notes-', notes, ...
                    '.mat');
    save(fname)
    fprintf('results saved to: \n %s \n', fname)
end
