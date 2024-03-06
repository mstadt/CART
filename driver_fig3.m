% Run simulation of clinical PKPD model
clear all;

%% set initial condition
CARTe_PB0 = 0; % CARTe in blood
CARTm_PB0 = 0; % CARTm in blood
CARTe_T0  = 1e-6; %0; %1e-100; % CARTe in tissue, seems to need to be nonzero to prevent numerical issues....
CARTm_T0  = 0; % CARTm in tissue
Cplx0     = 0; % CAR-Target Complexes
Tumor0    = 1e5; % tumor size
M0        = 12.1*2.5E9/0.117; % 
sBCMA0    = 0.5;

IC = [CARTe_PB0;CARTm_PB0;CARTe_T0;CARTm_T0;Cplx0;Tumor0;M0;sBCMA0];

%% set parameters
p = set_params('PKPD_clin');
% change params here
p.Tumor0 = Tumor0;

[params, ~] = pars2vector(p, 0);


%% time span
t0 = 0;
tf = 360; 
tspan = [t0,tf];

%% CART dose
doseCART_tot = 150e6; %50e6; %50e6; %50e6; % total number of cells in dose
dose_start = 0; % time to start dose
dose_time_hrs = 4; % time over dose is given (hrs)

%% ODE settings
options = odeset('RelTol',1.0e-12,'AbsTol',1e-12); % ode solver settings
%options = odeset('RelTol', 1.0e-12, 'AbsTol', 1e-12); % ode solver

%% Simulation
fprintf('vehicle simulation \n')
[t_veh,y_veh] = ode15s(@(t,y) modeqns_PKPD(t,y,params,...
                                    'doseCART', 0),...
                                    tspan, IC, options);
fprintf('treatment simulation \n')
[t_treat,y_treat] = run_dose_sim(doseCART_tot, dose_start,dose_time_hrs,...
                                    params, tspan, IC, options);
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

%% variables
figure(1)
nr = 4;
nc = 2;
clf;
subplot(nr,nc,1)
hold on
plot(t_veh,y_veh(:,1),'linewidth',lw,'color',c1)
plot(t_treat,y_treat(:,1),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTe_{PB}')
ylim([0, max(ymax, max(y_veh(:,1)))])
legend(labs)
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,2)
hold on
plot(t_veh,y_veh(:,2),'linewidth',lw,'color',c1)
plot(t_treat,y_treat(:,2),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTm_{PB}')
ylim([0, max(ymax, max(y_veh(:,2)))])
legend(labs)
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,3)
hold on
plot(t_veh,y_veh(:,3),'linewidth',lw,'color',c1)
plot(t_treat,y_treat(:,3),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTe_{T}')
ylim([0, max(ymax, max(y_veh(:,3)))])
legend(labs)
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,4)
hold on
plot(t_veh,y_veh(:,4),'linewidth',lw,'color',c1)
plot(t_treat,y_treat(:,4),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTm_{T}')
ylim([0, max(ymax, max(y_veh(:,4)))])
legend(labs)
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,5)
hold on
plot(t_veh,y_veh(:,5),'linewidth',lw,'color',c1)
plot(t_treat,y_treat(:,5),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('Cplx')
ylim([0, max(ymax, max(y_veh(:,5)))])
legend(labs)
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,6)
hold on
plot(t_veh,y_veh(:,6),'linewidth',lw,'color',c1)
plot(t_treat,y_treat(:,6),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('Tumor')
legend(labs)
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,7)
hold on
plot(t_veh,(y_veh(:,7)-M0) * 100/ M0,'linewidth',lw,'color',c1)
plot(t_treat,(y_treat(:,7)-M0) * 100/M0, 'linewidth',lw,'color',c2)
xlabel('t')
ylabel('M change')
legend(labs)
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,8)
hold on
plot(t_veh,(y_veh(:,8) - sBCMA0)*100/sBCMA0,'linewidth',lw,'color',c1)
plot(t_treat,(y_treat(:,8)-sBCMA0)*100/sBCMA0, 'linewidth',lw,'color',c2)
xlabel('t')
ylabel('sBMCA change')
legend(labs)
set(gca,'fontsize',f.gca)
grid on

%% treatment only
figure(2)
nr = 4;
nc = 2;
clf;
subplot(nr,nc,1)
plot(t_treat,y_treat(:,1),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTe_{PB}')
ylim([0, max(ymax, max(y_treat(:,1)))])
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,2)
plot(t_treat,y_treat(:,2),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTm_{PB}')
ylim([0, max(ymax, max(y_treat(:,2)))])
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,3)
plot(t_treat,y_treat(:,3),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTe_{T}')
ylim([0, max(ymax, max(y_treat(:,3)))])
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,4)
plot(t_treat,y_treat(:,4),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTm_{T}')
ylim([0, max(ymax, max(y_treat(:,4)))])
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,5)
plot(t_treat,y_treat(:,5),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('Cplx')
ylim([0, max(ymax, max(y_treat(:,5)))])
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,6)
plot(t_treat,y_treat(:,6),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('Tumor')
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,7)
plot(t_treat,(y_treat(:,7)-M0) * 100/M0, 'linewidth',lw,'color',c2)
xlabel('t')
ylabel('M change')
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,8)
plot(t_treat,(y_treat(:,8)-sBCMA0)*100/sBCMA0, 'linewidth',lw,'color',c2)
xlabel('t')
ylabel('sBMCA change')
set(gca,'fontsize',f.gca)
grid on

temp = sprintf('Treatment dose %d cells', doseCART_tot);
sgtitle(temp)


%% biomarkers
figure(3)
nr = 2;
nc = 2;
clf;
subplot(nr,nc,1)
plot(t_treat,y_treat(:,7), 'linewidth',lw,'color',c2)
xlabel('t')
ylabel('M')
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,2)
plot(t_treat,y_treat(:,8), 'linewidth',lw,'color',c2)
xlabel('t')
ylabel('sBCMA')
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,3)
plot(t_treat,(y_treat(:,7)-M0) * 100/M0, 'linewidth',lw,'color',c2)
xlabel('t')
ylabel('Change from baseline of serum M-protein (%)')
set(gca,'fontsize',f.gca)
grid on

subplot(nr,nc,4)
plot(t_treat,(y_treat(:,8)-sBCMA0)*100/sBCMA0, 'linewidth',lw,'color',c2)
xlabel('t')
ylabel('Change from baseline of serum BCMA (%)')
set(gca,'fontsize',f.gca)
grid on

temp = sprintf('Treatment dose %d cells', doseCART_tot);
sgtitle(temp)

%% tumor
figure(4)
clf;
plot(t_treat,y_treat(:,6)/Tumor0,'linewidth',lw,'color',c2)
xlabel('t')
ylabel('Tumor/Tumor_0')
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
