% Run simulation for Fig 3
clear all;

%% Dose
doseCART_tot = 50e6;

%% Set initial condition
CARTe_PB0 = 0; % CARTe in blood
CARTm_PB0 = 0; % CARTm in blood
CARTe_T0  = 0; %1e-6; %0; %1e-100; % CARTe in tissue, seems to need to be nonzero to prevent numerical issues....
CARTm_T0  = 0; % CARTm in tissue
Cplx0     = 0; % CAR-Target Complexes
Tumor0    = 1e5; % tumor size

IC_veh = [CARTe_PB0;CARTm_PB0;CARTe_T0;CARTm_T0;Cplx0;Tumor0];
IC_treat = [CARTe_PB0 + doseCART_tot;...
    CARTm_PB0;CARTe_T0;CARTm_T0;Cplx0;Tumor0];

%% set parameters
p = set_params('PKPD_clin');
% change params here
p.Rm = 1;

[params, ~] = pars2vector(p, 0);

%% time span
t0 = 0;
tf = 120;
tspan = [t0,tf];

%% ODE settings
options = odeset('RelTol',1e-9, 'AbsTol',1e-9); % ODE solver settings

%% Simulation
fprintf('start treatment simulation \n')
[t_treat,y_treat] = ode15s(@(t,y) modeqns_PKPD(t,y,params),...
                                    tspan, IC_treat, options);

% set small values to 0
small_vals = find(abs(y_treat)<0.01);
y_treat(small_vals) = 0; % set to 0

% % Negative values
% neg_vals = find(y_treat < 0);
% if max(abs(y_treat(neg_vals))) > 0.1
%     fprintf('warning: values went too negative!')
% else
%     y_treat = max(y_treat,0); % set all negative values to 0
% end


%% Plot results
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
xmax = tf;
%% Blood transgene
TransC=0.002; % from model code
CARTe_PB = y_treat(:,1);
CARTm_PB = y_treat(:,2);
FinalCARTPB = TransC * (CARTe_PB + CARTm_PB);
figure(1)
clf
plot(t_treat, FinalCARTPB, 'linewidth',lw,'color',c2)
xlabel('t (days)')
ylabel('Transgene copies/\mu g genomic DNA')
xlim([0,xmax])
ytick_values = [0,1e1, 1e2, 1e3, 1e4,1e5,1e6,1e7];
yticks(ytick_values);
set(gca,'Yscale','log','fontsize',f.gca)
grid on

%% CART in each compartment
figure(2)
nr = 2; nc = 2;
subplot(nr,nc,1)
plot(t_treat, CARTe_PB, 'linewidth',lw,'color',c2)
xlabel('t (days)')
ylabel('CARTePB')
xlim([0,xmax])
%ytick_values = [0,1e1, 1e2, 1e3, 1e4,1e5,1e6,1e7];
%yticks(ytick_values);
set(gca,'Yscale','log','fontsize',f.gca)
grid on

subplot(nr,nc,2)
plot(t_treat, CARTm_PB, 'linewidth',lw,'color',c2)
xlabel('t (days)')
ylabel('CARTmPB')
xlim([0,xmax])
%ytick_values = [0,1e1, 1e2, 1e3, 1e4,1e5,1e6,1e7];
%yticks(ytick_values);
set(gca,'Yscale','log','fontsize',f.gca)
grid on

subplot(nr,nc,3)
plot(t_treat, y_treat(:,3), 'linewidth',lw,'color',c2)
xlabel('t (days)')
ylabel('CARTeT')
xlim([0,xmax])
%ytick_values = [0,1e1, 1e2, 1e3, 1e4,1e5,1e6,1e7];
%yticks(ytick_values);
set(gca,'Yscale','log','fontsize',f.gca)
grid on

subplot(nr,nc,4)
plot(t_treat, y_treat(:,4), 'linewidth',lw,'color',c2)
xlabel('t (days)')
ylabel('CARTmT')
xlim([0,xmax])
%ytick_values = [0,1e1, 1e2, 1e3, 1e4,1e5,1e6,1e7];
%yticks(ytick_values);
set(gca,'Yscale','log','fontsize',f.gca)
grid on

