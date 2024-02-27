% Run simulation of PD model
clear all;

%% set parameters
p = set_params('PD');
[params, ~] = pars2vector(p, 0);

%% set initial condition
Cmplx0 = 0;
Nt0    = 1e5;
Ne0    = 0;

IC = [Cmplx0; Nt0; Ne0];

%% time span
t0 = 0;
tf = 500;
tspan = [t0,tf];

%% ODE options
options = odeset('RelTol',1.0e-6,'AbsTol',1e-9); % ode solver settings

%% Simulation
fprintf('solving ODEs \n')
[t,y] = ode15s(@(t,y) modeqns_PD(t,y,params),...
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

% variables
figure(1)
clf;
subplot(1,3,1)
plot(t,y(:,1),'linewidth',lw,'color',c1)
xlabel('t')
ylabel('Cmplx')
set(gca,'fontsize',f.gca)
grid on

subplot(1,3,2)
plot(t,y(:,2),'linewidth',lw,'color',c1)
xlabel('t')
ylabel('Nt')
set(gca,'fontsize',f.gca)
grid on

subplot(1,3,3)
plot(t,y(:,3),'linewidth',lw,'color',c1)
xlabel('t')
ylabel('Ne')
set(gca,'fontsize',f.gca)
grid on

%% save simulations options
save_sim = input('Do you want to save the simulation? (0 - no/1 - yes) ');
if save_sim
    notes = input('notes: ');
    fname = strcat('./simPD/', date, '_driver_PD', ...
                    '_notes-', notes, ...
                    '.mat');
    save(fname)
    fprintf('results saved to: \n %s \n', fname)
end