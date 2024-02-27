% Run simulation of clinical PKPD model
clear all;

%% set parameters
p = set_params('PKPD_clin');
[params, ~] = pars2vector(p, 0);

%% set initial condition
CARTe_PB0 = 0; % CARTe in blood
CARTm_PB0 = 0; % CARTm in blood
CARTe_T0  = 1e-100; % CARTe in tissue, seems to need to be nonzero to prevent numerical issues....
CARTm_T0  = 0; % CARTm in tissue
Cplx0     = 0; % CAR-Target Complexes
Tumor0    = 1e5; % tumor size

IC = [CARTe_PB0;CARTm_PB0;CARTe_T0;CARTm_T0;Cplx0;Tumor0];

%% time span
t0 = 0;
tf = 50;
tspan = [t0,tf];

%% CART dose
doseCART = 0; 

%% ODE settings
options = odeset('RelTol',1.0e-12,'AbsTol',1e-12); % ode solver settings

%% Simulation
fprintf('vehicle simulation \n')
% vehicle simulation
[t,y] = ode15s(@(t,y) modeqns_PKPD(t,y,params,...
                                    'doseCART', 0),...
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

%ymin = 0;
ymax = 1;

% variables
figure(1)
nr = 3;
nc = 2;
clf;
subplot(nr,nc,1)
plot(t,y(:,1),'linewidth',lw,'color',c1)
xlabel('t')
ylabel('CARTe_{PB}')
ylim([0, max(ymax, max(y(:,1)))])
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,2)
plot(t,y(:,2),'linewidth',lw,'color',c1)
xlabel('t')
ylabel('CARTm_{PB}')
ylim([0, max(ymax, max(y(:,2)))])
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,3)
plot(t,y(:,3),'linewidth',lw,'color',c1)
xlabel('t')
ylabel('CARTe_{T}')
ylim([0, max(ymax, max(y(:,3)))])
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,4)
plot(t,y(:,4),'linewidth',lw,'color',c1)
xlabel('t')
ylabel('CARTm_{T}')
ylim([0, max(ymax, max(y(:,4)))])
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,5)
plot(t,y(:,5),'linewidth',lw,'color',c1)
xlabel('t')
ylabel('Cplx')
ylim([0, max(ymax, max(y(:,5)))])
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,6)
plot(t,y(:,6),'linewidth',lw,'color',c1)
xlabel('t')
ylabel('Tumor')
set(gca,'fontsize',f.gca)
grid on

%% tumor volume and CART cells in blood
% get extracted data from Ruiz-Martinez Fig 2b
dat = load('./data/fig2b_data.mat');
figure(2);
ms = 15;
cmap2 = parula(6);
c1 = cmap2(2,:);
c2 = cmap2(5,:);
clf;
nr = 1; nc = 2;
% Plot data
subplot(nr,nc,1)
hold on
plot(dat.datTV_veh(:,1), dat.datTV_veh(:,2),...
                'linestyle','none',...
                'marker', 'o', 'markersize', ms,...
                'color', c1,'markerfacecolor', c1)
plot(dat.datTV_treat(:,1), dat.datTV_treat(:,2),...
                'linestyle','none',...
                'marker', 'o', 'markersize', ms,...
                'color', c2,'markerfacecolor', c2)
xlabel('Time (day)')
ylabel('Tumor volume (mm^3)')
set(gca,'fontsize',f.gca)
grid on
subplot(nr,nc,2)
hold on
plot(dat.datCART_veh(:,1), dat.datCART_veh(:,2),...
                'linestyle','none',...
                'marker', 'o', 'markersize', ms,...
                'color', c1,'markerfacecolor', c1)
plot(dat.datCART_treat(:,1), dat.datCART_treat(:,2),...
                'linestyle','none',...
                'marker', 'o', 'markersize', ms,...
                'color', c2,'markerfacecolor', c2)
xlabel('Time (day)')
ylabel('CAR-T Cells in Blood (#/\muL)')
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