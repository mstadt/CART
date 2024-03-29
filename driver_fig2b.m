% Run simulation of preclinical PKPD model with and without drugs
clear all;

%% set parameters
p = set_params('PKPD_preclin');
% change params here (option)
p.Kg_tumor = 0.105;

% density_CAR
p.Ag_CAR = 20000; % option in code

p.Rm = 2e-2; %2e-2; %2e-2; %2e-2; 



[params, ~] = pars2vector(p, 0);

%% set initial condition
CARTe_PB0 = 0; % CARTe in blood
CARTm_PB0 = 0; % CARTm in blood
CARTe_T0  = 0; % CARTe in tissue, seems to need to be nonzero to prevent numerical issues....
CARTm_T0  = 0; % CARTm in tissue
Cplx0     = 0; % CAR-Target Complexes
Tumor0    = 2.5e9; %6e7; %1e5; % tumor size

IC = [CARTe_PB0;CARTm_PB0;CARTe_T0;CARTm_T0;Cplx0;Tumor0];

%% time span
t0 = 0;
tf_veh = 35;
tf_treat = 70;
tspan_veh = [t0,tf_veh];
tspan_treat = [t0;tf_treat];

%% CART dose
doseCART_tot = 5e6; % total number of cells in dose


%% ODE settings
options = odeset('RelTol',1.0e-12,'AbsTol',1e-16); % ode solver settings

%% Simulation
fprintf('vehicle simulation \n')
% vehicle simulation
[t_veh, y_veh] = ode15s(@(t,y) modeqns_PKPD(t,y,params),...
                                    tspan_veh, IC, options);
fprintf('treatment simulation \n')
IC(1) = doseCART_tot/p.Vb; % add dose/L
[t_treat, y_treat] = ode15s(@(t,y) modeqns_PKPD(t,y,params),...
                                    tspan_treat, IC, options);
% [t_treat,y_treat] = run_dose_sim(doseCART_tot, dose_start,dose_time_hrs,...
%                                     params, tspan_treat, IC, options);
fprintf('sim finished \n')

%% Make figures
fprintf('making figs \n')
% figure specs
lw = 4;
f.xlab = 16; f.ylab = 16; f.title = 18;
f.leg = 16; f.gca = 12;
cmap = parula(6);
c1 = cmap(1,:);
c2 = cmap(4,:);
ls1 = '-'; ls2 = '-';
cgraymap = gray(5);
cgray = cgraymap(3,:);
lwgray = 2; lsgray = '--';

labs = {'vehicle', 'treatment'};

% variables
% Vehicle
figure(1)
nr = 3;
nc = 2;
clf;
subplot(nr,nc,1)
hold on
plot(t_veh,y_veh(:,1),'linewidth',lw,'color',c1)
plot(t_treat,max(1e-5,y_treat(:,1)),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTe_{PB}')
legend(labs)
ylim([10^0, 10^8])
yticks([10^0, 10^2, 10^4, 10^6, 10^8])
set(gca,'fontsize',f.gca,'Yscale','log')
grid on
hold off

subplot(nr,nc,2)
hold on
plot(t_veh,y_veh(:,2),'linewidth',lw,'color',c1)
plot(t_treat,y_treat(:,2),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTm_{PB}')
legend(labs)
%ymax = max([1;y_veh(:,2);y_treat(:,2)]);
ylim([10^0, 10^8])
yticks([10^0, 10^2, 10^4, 10^6, 10^8])
set(gca,'fontsize',f.gca,'Yscale','log')
grid on
hold off

subplot(nr,nc,3)
hold on
plot(t_veh,y_veh(:,3),'linewidth',lw,'color',c1)
plot(t_treat,y_treat(:,3),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTe_{T}')
legend(labs)
ylim([10^0, 10^8])
yticks([10^0, 10^2, 10^4, 10^6, 10^8])
set(gca,'fontsize',f.gca,'Yscale','log')
grid on
hold off

subplot(nr,nc,4)
hold on
plot(t_veh,y_veh(:,4),'linewidth',lw,'color',c1)
plot(t_treat,y_treat(:,4),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTm_{T}')
legend(labs)
%ymax = max([1;y_veh(:,4);y_treat(:,4)]);
ylim([10^0, 10^8])
yticks([10^0, 10^2, 10^4, 10^6, 10^8])
set(gca,'fontsize',f.gca,'Yscale','log')
grid on
hold off

subplot(nr,nc,5)
hold on
plot(t_veh,y_veh(:,5),'linewidth',lw,'color',c1)
plot(t_treat,y_treat(:,5),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('Cplx')
legend(labs)
%ymax = max([1;y_veh(:,5);y_treat(:,5)]);
%ylim([0, ymax])
set(gca,'fontsize',f.gca)
grid on
hold off

subplot(nr,nc,6)
hold on
plot(t_veh,y_veh(:,6),'linewidth',lw,'color',c1)
plot(t_treat,y_treat(:,6),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('Tumor')
legend(labs)
set(gca,'fontsize',f.gca)
grid on
hold off



%% tumor volume and CART cells in blood
% get extracted data from Ruiz-Martinez Fig 2b
dat = load('./data/fig2b_data.mat');
figure(2);
ms = 15;
clf;
nr = 1; nc = 2;
% Tumor volume
subplot(nr,nc,1)
hold on
% vehicle
plot(dat.datTV_veh(:,1), dat.datTV_veh(:,2),...
                'linestyle','none',...
                'marker', 'o', 'markersize', ms,...
                'color', c1,'markerfacecolor', c1, ...
                'HandleVisibility', 'off')

TC2Vol =100/Tumor0; % cells 2 volume conversion
plot(t_veh,y_veh(:,6)*TC2Vol,'linewidth',lw,'color',c1)
% treatment
plot(dat.datTV_treat(:,1), dat.datTV_treat(:,2),...
                'linestyle','none',...
                'marker', 'o', 'markersize', ms,...
                'color', c2,'markerfacecolor', c2, ...
                'HandleVisibility', 'off')
plot(t_treat,y_treat(:,6)*TC2Vol,'linewidth',lw,'color',c2)
legend(labs,'fontsize',18)
xlabel('Time (day)')
ylabel('Tumor volume (mm^3)')
ylim([0,4000])
set(gca,'fontsize',f.gca)
grid on
hold off

% CART in blood
subplot(nr,nc,2)
hold on
% vehicle
plot(dat.datCART_veh(:,1), dat.datCART_veh(:,2),...
                'linestyle','none',...
                'marker', 'o', 'markersize', ms,...
                'color', c1,'markerfacecolor', c1, ...
                'HandleVisibility', 'off')

CART_PB_tot = (y_veh(:,1) + y_veh(:,2));%(p.Vb * 1000); % # cell/muL
plot(t_veh,CART_PB_tot,'linewidth',lw,'color',c1)
% treatment
plot(dat.datCART_treat(:,1), dat.datCART_treat(:,2),...
                'linestyle','none',...
                'marker', 'o', 'markersize', ms,...
                'color', c2,'markerfacecolor', c2, ...
                'HandleVisibility', 'off')
cf = 1e-3;
CART_PB_tot = (y_treat(:,1) + y_treat(:,2))*cf; %(p.Vb * 1000); % # cell/muL
plot(t_treat,CART_PB_tot,'linewidth',lw,'color',c2)
xlabel('Time (day)')
ylabel('CAR-T Cells in Blood (#/\muL)')
set(gca,'fontsize',f.gca)
xlim([0,40])
ylim([0,1000])
legend(labs,'fontsize',18)
grid on
hold off

AddLetters2Plots(figure(2),{'A','B'},'HShift', -0.06, 'VShift', -0.06, ...
                'fontsize', 22)

%% Plot dosing
figure(3)
lw = 4;
nr = 2; nc = 2;
subplot(nr,nc,1)
plot(t_treat, y_treat(:, 1) * p.Vb,'linewidth',lw, 'color',c2)
xlabel('t (day)')
ylabel('CARTe_{PB} (number of cells)')
ylim([10^0,10^11])
set(gca,'fontsize',f.gca, 'Yscale','log')
grid on

subplot(nr,nc,2)
plot(t_treat, y_treat(:, 3) * p.Vt,'linewidth',lw, 'color',c2)
xlabel('t (day)')
ylabel('CARTe_{T} (number of cells)')
ylim([10^0,10^11])
set(gca,'fontsize',f.gca, 'YScale','log')
grid on

subplot(nr,nc,3)
plot(t_treat, y_treat(:, 2) * p.Vb,'linewidth',lw, 'color',c2)
xlabel('t (day)')
ylabel('CARTm_{PB} (number of cells)')
ylim([10^0,10^11])
set(gca,'fontsize',f.gca, 'Yscale','log')
grid on

subplot(nr,nc,4)
plot(t_treat, y_treat(:, 4) * p.Vt,'linewidth',lw, 'color',c2)
xlabel('t (day)')
ylabel('CARTm_{T} (number of cells)')
ylim([10^0,10^11])
set(gca,'fontsize',f.gca, 'YScale','log')
grid on

AddLetters2Plots(figure(3),{'A','B','C','D'},'HShift', -0.06, 'VShift', -0.06, ...
                'fontsize', 22)


%% Cplx
figure(4)
plot(t_treat*24, y_treat(:,5), 'linewidth', 3, 'color', c2)
xlabel('t (day)')
ylabel('Cplx')
set(gca,'fontsize',f.gca)
grid on


y=y_treat;
t = t_treat;
run('compute_alg_eqns.m')
