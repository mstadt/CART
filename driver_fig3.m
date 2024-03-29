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
Tumor0    = 2.5e9; %1e5; % tumor size

IC = [CARTe_PB0;CARTm_PB0;CARTe_T0;CARTm_T0;Cplx0;Tumor0];

%% set parameters
p = set_params('PKPD_clin');
% change params here
p.Rm = 2e-3;

[params, ~] = pars2vector(p, 0);


%% time span
t0 = 0;
tf = 360; 
tspan = [t0,tf];

%% CART dose
doseCART_tot = 800e6; %500e6; %; % total number of cells in dose
%dose_start = 0; % time to start dose
%dose_time_hrs = 4; % time over dose is given (hrs)

%% ODE settings
options = odeset('RelTol',1.0e-12,'AbsTol',1e-12); % ode solver settings

%% Simulation
% fprintf('vehicle simulation \n')
% [t_veh,y_veh] = ode15s(@(t,y) modeqns_PKPD(t,y,params),...
%                                     tspan, IC, options);
fprintf('treatment simulation \n')
IC(1) = doseCART_tot;
[t_treat, y_treat] = ode15s(@(t,y) modeqns_PKPD(t,y,params),...
                                    tspan, IC, options);
fprintf('sim finished \n')

if do_biomarkers
    M0        = 1; %12.1*2.5E9/0.117; % 
    sBCMA0    = 1; %0.5;
    % biomarkers simulation
    fprintf('get biomarkers \n')
    IC_bm = [M0;sBCMA0];
    tTumor = t_treat;
    Tumor = y_treat(:,6);
    % biomarker sim parameters
    p_bm = set_params('biomarkers');
    p_bm.Tumor0 = Tumor(1);
    [pars_bm,~] = pars2vector(p_bm,0);
    
    % Check for uniqueness
    % Find the unique values and their counts
    [uniqueValues, ~, index] = unique(tTumor);
    counts = histcounts(tTumor, uniqueValues);
    
    % Find the values that appear more than once
    nonUniqueValues = uniqueValues(counts > 1);
    
    if ~isempty(nonUniqueValues)
        for ii = 1:length(nonUniqueValues)
            val = nonUniqueValues(ii);
            inds = find(tTumor == val);
            % remove values after first value
            for jj = 2:length(inds)
                id = inds(jj);
                tTumor(id) = [];
                Tumor(id) = [];
            end
        end
    end
    
    
    [t_bm,y_bm] = ode15s(@(t,y) modeqns_biomarkers(t,y,...
                                        interp1(tTumor, Tumor, t),...
                                        pars_bm),...
                                        tspan, IC_bm, options);
end



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

% %% variables
% figure(1)
% nr = 3;
% nc = 2;
% clf;
% subplot(nr,nc,1)
% hold on
% plot(t_veh,y_veh(:,1),'linewidth',lw,'color',c1)
% plot(t_treat,y_treat(:,1),'linewidth',lw,'color',c2)
% xlabel('t')
% ylabel('CARTe_{PB}')
% ylim([0, max(ymax, max(y_veh(:,1)))])
% legend(labs)
% set(gca,'fontsize',f.gca)
% grid on
% 
% subplot(nr,nc,2)
% hold on
% plot(t_veh,y_veh(:,2),'linewidth',lw,'color',c1)
% plot(t_treat,y_treat(:,2),'linewidth',lw,'color',c2)
% xlabel('t')
% ylabel('CARTm_{PB}')
% ylim([0, max(ymax, max(y_veh(:,2)))])
% legend(labs)
% set(gca,'fontsize',f.gca)
% grid on
% 
% subplot(nr,nc,3)
% hold on
% plot(t_veh,y_veh(:,3),'linewidth',lw,'color',c1)
% plot(t_treat,y_treat(:,3),'linewidth',lw,'color',c2)
% xlabel('t')
% ylabel('CARTe_{T}')
% ylim([0, max(ymax, max(y_veh(:,3)))])
% legend(labs)
% set(gca,'fontsize',f.gca)
% grid on
% 
% subplot(nr,nc,4)
% hold on
% plot(t_veh,y_veh(:,4),'linewidth',lw,'color',c1)
% plot(t_treat,y_treat(:,4),'linewidth',lw,'color',c2)
% xlabel('t')
% ylabel('CARTm_{T}')
% ylim([0, max(ymax, max(y_veh(:,4)))])
% legend(labs)
% set(gca,'fontsize',f.gca)
% grid on
% 
% subplot(nr,nc,5)
% hold on
% plot(t_veh,y_veh(:,5),'linewidth',lw,'color',c1)
% plot(t_treat,y_treat(:,5),'linewidth',lw,'color',c2)
% xlabel('t')
% ylabel('Cplx')
% ylim([0, max(ymax, max(y_veh(:,5)))])
% legend(labs)
% set(gca,'fontsize',f.gca)
% grid on
% 
% subplot(nr,nc,6)
% hold on
% plot(t_veh,y_veh(:,6),'linewidth',lw,'color',c1)
% plot(t_treat,y_treat(:,6),'linewidth',lw,'color',c2)
% xlabel('t')
% ylabel('Tumor')
% legend(labs)
% set(gca,'fontsize',f.gca)
% grid on


%% treatment only
figure(2)
nr = 3;
nc = 2;
clf;
subplot(nr,nc,1)
plot(t_treat,max(1e-16,y_treat(:,1)),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTe_{PB}')
ylim([0, max(ymax, max(y_treat(:,1)))])
set(gca,'fontsize',f.gca,'YScale','log')
grid on

subplot(nr,nc,2)
plot(t_treat,y_treat(:,2),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTm_{PB}')
ylim([0, max(ymax, max(y_treat(:,2)))])
set(gca,'fontsize',f.gca,'YScale','log')
grid on

subplot(nr,nc,3)
plot(t_treat,y_treat(:,3),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTe_{T}')
ylim([0, max(ymax, max(y_treat(:,3)))])
set(gca,'fontsize',f.gca,'YScale','log')
grid on

subplot(nr,nc,4)
plot(t_treat,y_treat(:,4),'linewidth',lw,'color',c2)
xlabel('t')
ylabel('CARTm_{T}')
ylim([0, max(ymax, max(y_treat(:,4)))])
set(gca,'fontsize',f.gca,'YScale','log')
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

temp = sprintf('Treatment dose %d cells', doseCART_tot);
sgtitle(temp)


if do_biomarkers
    %% biomarkers
    figure(3)
    nr = 2;
    nc = 2;
    clf;
    subplot(nr,nc,1)
    plot(t_bm,y_bm(:,1), 'linewidth',lw,'color',c2)
    xlabel('t')
    ylabel('M')
    set(gca,'fontsize',f.gca)
    grid on
    
    subplot(nr,nc,2)
    plot(t_bm,y_bm(:,2), 'linewidth',lw,'color',c2)
    xlabel('t')
    ylabel('sBCMA')
    set(gca,'fontsize',f.gca)
    grid on
    
%     subplot(nr,nc,3)
%     plot(t_bm,(y_bm(:,1)-M0) * 100/M0, 'linewidth',lw,'color',c2)
%     xlabel('t')
%     ylabel('Change from baseline of serum M-protein (%)')
%     set(gca,'fontsize',f.gca)
%     grid on
%     
%     subplot(nr,nc,4)
%     plot(t_bm,(y_bm(:,2)-sBCMA0)*100/sBCMA0, 'linewidth',lw,'color',c2)
%     xlabel('t')
%     ylabel('Change from baseline of serum BCMA (%)')
%     set(gca,'fontsize',f.gca)
%     grid on
    
    temp = sprintf('Treatment dose %d cells', doseCART_tot);
    sgtitle(temp)
end

%% tumor
figure(4)
clf;
hold on
%plot(t_veh,y_veh(:,6)/Tumor0,'linewidth',lw,'color',c1)
plot(t_treat,y_treat(:,6)/Tumor0,'linewidth',lw,'color',c2)
xlabel('t')
ylabel('Tumor/Tumor_0')
%legend('Vehicle', 'Treatment')
set(gca,'fontsize',f.gca)
grid on

%% Blood transgene
TransC=0.002; % from model code
CARTe_PB = y_treat(:,1);
CARTm_PB = y_treat(:,2);
FinalCARTPB = TransC * (CARTe_PB + CARTm_PB);
figure(5)
%clf
plot(t_treat, FinalCARTPB, 'linewidth',lw)
%plot(t_treat, FinalCARTPB, 'linewidth',lw,'color',c2)
xlabel('t (days)')
ylabel('Transgene copies/\mu g genomic DNA')
xlim([0,70])
ylim([10^1, 10^7])
set(gca,'fontsize',f.gca, 'Yscale', 'log')
grid on

%% Cplx
figure(6)
clf;
plot(t_treat, y_treat(:,5),'linewidth',lw,'color',c2)
xlabel('t (days)')
ylabel('Cplx')
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

y=y_treat;
t=t_treat;

