% Run simulation of preclinical PKPD model with and without drugs
clear all;

%% set parameters
p = set_params('PKPD_preclin');
pclin = set_params('PKPD_clin');
% change params here (option)
p = pclin;
p.Kg_tumor = 0.1; %0.2; %0.088;
% p.K12 = pclin.K12;
% p.Kexp_max = pclin.Kexp_max;
% p.EC50_exp = pclin.EC50_exp;


[params, ~] = pars2vector(p, 0);

%% set initial condition
CARTe_PB0 = 0; % CARTe in blood
CARTm_PB0 = 0; % CARTm in blood
CARTe_T0  = 0; % CARTe in tissue, seems to need to be nonzero to prevent numerical issues....
CARTm_T0  = 0; % CARTm in tissue
Cplx0     = 0; % CAR-Target Complexes
Tumor0    = 1e5; % tumor size

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
IC(1) = doseCART_tot; % add dose
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
f.leg = 16; f.gca = 18;
cmap = parula(6);
c1 = cmap(2,:);
c2 = cmap(5,:);
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
plot(t_treat,y_treat(:,1),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTe_{PB}')
legend(labs)
ymax = max([1;y_veh(:,1);y_treat(:,1)]);
ylim([0, ymax])
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
ymax = max([1;y_veh(:,2);y_treat(:,2)]);
ylim([0, ymax])
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
ymax = max([1;y_veh(:,3);y_treat(:,3)]);
ylim([0, ymax])
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
ymax = max([1;y_veh(:,4);y_treat(:,4)]);
ylim([0, ymax])
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
ymax = max([1;y_veh(:,5);y_treat(:,5)]);
ylim([0, ymax])
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

TC2Vol = 100/1e5; % cells 2 volume conversion
plot(t_veh,y_veh(:,6)*TC2Vol,'linewidth',lw,'color',c1)
% treatment
plot(dat.datTV_treat(:,1), dat.datTV_treat(:,2),...
                'linestyle','none',...
                'marker', 'o', 'markersize', ms,...
                'color', c2,'markerfacecolor', c2, ...
                'HandleVisibility', 'off')
plot(t_treat,y_treat(:,6)*TC2Vol,'linewidth',lw,'color',c2)
legend(labs)
xlabel('Time (day)')
ylabel('Tumor volume (mm^3)')
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

CART_PB_tot = (y_veh(:,1) + y_veh(:,2))/5000;%(p.Vb * 1000); % # cell/muL
plot(t_veh,CART_PB_tot,'linewidth',lw,'color',c1)
% treatment
plot(dat.datCART_treat(:,1), dat.datCART_treat(:,2),...
                'linestyle','none',...
                'marker', 'o', 'markersize', ms,...
                'color', c2,'markerfacecolor', c2, ...
                'HandleVisibility', 'off')
CART_PB_tot = (y_treat(:,1) + y_treat(:,2))/5000; %(p.Vb * 1000); % # cell/muL
plot(t_treat,CART_PB_tot,'linewidth',lw,'color',c2)
legend(labs)
xlabel('Time (day)')
ylabel('CAR-T Cells in Blood (#/\muL)')
set(gca,'fontsize',f.gca)
xlim([0,40])
grid on
hold off

%% Plot dosing
figure(3)
ind = find(t_treat > 240/24, 1, 'first');

subplot(1,2,1)
plot(t_treat(1:ind) * 24, y_treat(1:ind, 1) * p.Vb,'linewidth',3, 'color',c2)
xlabel('t (hrs)')
ylabel('CARTe_{PB} (number of cells)')
set(gca,'fontsize',f.gca)
grid on

subplot(1,2,2)
plot(t_treat(1:ind) * 24, y_treat(1:ind, 3) * p.Vt,'linewidth',3, 'color',c2)
xlabel('t (hrs)')
ylabel('CARTe_{T} (number of cells)')
set(gca,'fontsize',f.gca)
grid on



%% save simulations options
% save_sim = input('Do you want to save the simulation? (0 - no/1 - yes) ');
% if save_sim
%     notes = input('notes: ');
%     fname = strcat('./simPD/', date, '_driver_PD', ...
%                     '_notes-', notes, ...
%                     '.mat');
%     save(fname)
%     fprintf('results saved to: \n %s \n', fname)
% end

