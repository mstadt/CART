% Run simulation of preclinical PKPD model
clear all;

%% set parameters
p = set_params('PKPD_preclin');
% change params here (option)

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
tf = 100;
tspan = [t0,tf];

%% CART dose
doseCART_tot = 5e6; % total number of cells in dose
dose_start = 0;
dose_time_hrs = 4; % time in hours


%% ODE settings
options = odeset('RelTol',1.0e-12,'AbsTol',1e-16); % ode solver settings

%% Simulation
fprintf('vehicle simulation \n')
% vehicle simulation
[t_veh, y_veh] = ode15s(@(t,y) modeqns_PKPD(t,y,params,...
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
set(gca,'fontsize',f.gca)
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
set(gca,'fontsize',f.gca)
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
set(gca,'fontsize',f.gca)
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
set(gca,'fontsize',f.gca)
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

CART_PB_tot = (y_veh(:,1) + y_veh(:,2))/(p.Vb * 1000); % # cell/muL
plot(t_veh,CART_PB_tot,'linewidth',lw,'color',c1)
% treatment
plot(dat.datCART_treat(:,1), dat.datCART_treat(:,2),...
                'linestyle','none',...
                'marker', 'o', 'markersize', ms,...
                'color', c2,'markerfacecolor', c2, ...
                'HandleVisibility', 'off')
CART_PB_tot = y_treat(:,1)/(p.Vb * 1000) %(y_treat(:,1) + y_treat(:,2))/(p.Vb * 1000); % # cell/muL
plot(t_treat,CART_PB_tot,'linewidth',lw,'color',c2)
legend(labs)
xlabel('Time (day)')
ylabel('CAR-T Cells in Blood (#/\muL)')
set(gca,'fontsize',f.gca)
grid on
hold off



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


%--------------------------
% functions used
%--------------------------
function [t,y] = run_dose_sim(doseCART_tot, dose_start,dose_time_hrs, params, tspan, IC, options)
    doseCART = doseCART_tot/(dose_time_hrs/24);
    if dose_start > tspan(1)
        % start with no dose to start
        t0 = tspan(1);
        tf = dose_start;
        IC1 = IC;
        [t1,y1] = ode15s(@(t,y) modeqns_PKPD(t,y,params,...
                                            'doseCART', 0),...
                                            [t0,tf], IC1, options);
        % add dose
        IC2 = y1(end,:);
        t0 = dose_start;
        tf = min(dose_start + dose_time_hrs/24, tspan(2));
        [t2,y2] = ode15s(@(t,y) modeqns_PKPD(t,y,params,...
                                            'doseCART', doseCART),...
                                            [t0,tf], IC2, options);

        t = [t1;t2];
        y = [y1;y2];
        if tspan(2) > (dose_start + dose_time_hrs/24)
            IC3 = y2(end,:);
            t0 = t2(end);
            tf = tspan(2);
            [t3,y3] = ode15s(@(t,y) modeqns_PKPD(t,y,params,...
                                            'doseCART', 0),...
                                            [t0,tf], IC3, options);
            t = [t;t3];
            y = [y;y3];
        end
    elseif dose_start == tspan(1)
        % start with dose
        IC1 = IC;
        t0 = tspan(1);
        tf = dose_time_hrs/24;
        [t1,y1] = ode15s(@(t,y) modeqns_PKPD(t,y,params,...
                                            'doseCART', doseCART),...
                                            [t0,tf], IC1, options);
        t = t1;
        y = y1;
       if tspan(2) > (dose_start + dose_time_hrs/24)
            IC2 = y1(end,:);
            t0 = t1(end);
            tf = tspan(2);
            [t2,y2] = ode15s(@(t,y) modeqns_PKPD(t,y,params,...
                                            'doseCART', 0),...
                                            [t0,tf], IC2, options);
            t = [t;t2];
            y = [y;y2];
        end
    else
        fprintf('dose_start: %f, tspan(1): %f \n', dose_start, tspan(1))
        error('dose_start before simulation start')
    end
end