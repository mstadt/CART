clear all;

% plot results for clinical model
% results made from driver_fig3_fromcode.m
f1 = './simPKPD/20-Mar-2024_driver_fig3_fromcode_dose-50000000_notes-reduceKkill.mat'; './simPKPD/20-Mar-2024_driver_fig3_fromcode_dose-50000000_notes-newpars.mat';%'./simPKPD/20-Mar-2024_driver_fig3_fromcode_dose-50000000_notes-dose1.mat';
f1b = './simPKPD/20-Mar-2024_driver_fig3_fromcode_dose-50000000_notes-origKkill.mat';
f2 = './simPKPD/20-Mar-2024_driver_fig3_fromcode_dose-150000000_notes-otherpar.mat';%'./simPKPD/20-Mar-2024_driver_fig3_fromcode_dose-150000000_notes-dose2.mat';
f3 = './simPKPD/20-Mar-2024_driver_fig3_fromcode_dose-450000000_notes-otherpar.mat'; %'./simPKPD/20-Mar-2024_driver_fig3_fromcode_dose-450000000_notes-dose3.mat';
f4 = './simPKPD/20-Mar-2024_driver_fig3_fromcode_dose-800000000_notes-otherpar.mat';%'./simPKPD/20-Mar-2024_driver_fig3_fromcode_dose-800000000_notes-dose4.mat';

dat1 = load(f1);
dat1b = load(f1b);
dat2 = load(f2);
dat3 = load(f3);
dat4 = load(f4);

% labels
lab1 = '50e6 - reduced Kkill\_max';
lab1b = '50e6 - original Kkill\_max';
lab2 = '150e6';
lab3 = '450e6';
lab4 = '800e6';

%% plot results
lw = 5;
f.xlab = 16; f.ylab = 16; f.title = 18;
f.leg = 16; f.gca = 18;
cmap = parula(5);
ls1 = '-'; ls2 = '-';
cgraymap = gray(5);
cgray = cgraymap(3,:);
lwgray = 2; lsgray = '--';
labs = {lab1,lab1b,lab2,lab3,lab4};

%% Biomarkers
M0 = dat1.M0;
B0 = dat1.B0;

figure(1)
clf;
nr = 1; nc = 2;

% M
subplot(nr,nc,1)
hold on
% dat 1
M = dat1.y_treat(:,7);
Mchange = (M - M0)/M0 * 100;
plot(dat1.t_treat, Mchange, 'linewidth', lw,'color',cmap(1,:))
M = dat1b.y_treat(:,7);
Mchange = (M - M0)/M0 * 100;
plot(dat1b.t_treat, Mchange, 'linewidth', lw,'color',cmap(1,:),'linestyle',':')
M = dat2.y_treat(:,7);
Mchange = (M - M0)/M0 * 100;
plot(dat2.t_treat, Mchange, 'linewidth', lw,'color',cmap(2,:))
M = dat3.y_treat(:,7);
Mchange = (M - M0)/M0 * 100;
plot(dat3.t_treat, Mchange, 'linewidth', lw,'color',cmap(3,:))
M = dat4.y_treat(:,7);
Mchange = (M - M0)/M0 * 100;
plot(dat4.t_treat, Mchange, 'linewidth', lw,'color',cmap(4,:))
xlim([0,180])
%ylim([-40,15])
xlabel('t (days)')
ylabel('Change from baseline of serum M-protein')
legend(labs)
set(gca,'fontsize',f.gca)
grid on

% B
subplot(nr,nc,2)
hold on
% dat 1
B = dat1.y_treat(:,8);
Bchange = (B - B0)/B0 * 100;
plot(dat1.t_treat, Bchange, 'linewidth', lw,'color',cmap(1,:))
B = dat1b.y_treat(:,8);
Bchange = (B - B0)/B0 * 100;
plot(dat1b.t_treat, Bchange, 'linewidth', lw,'color',cmap(1,:),'linestyle',':')
B = dat2.y_treat(:,8);
Bchange = (B - B0)/B0 * 100;
plot(dat2.t_treat, Bchange, 'linewidth', lw,'color',cmap(2,:))
B = dat3.y_treat(:,8);
Bchange = (B - B0)/B0 * 100;
plot(dat3.t_treat, Bchange, 'linewidth', lw,'color',cmap(3,:))
B = dat4.y_treat(:,8);
Bchange = (B - B0)/B0 * 100;
plot(dat4.t_treat, Bchange, 'linewidth', lw,'color',cmap(4,:))
xlim([0,180])
%ylim([-40,15])
xlabel('t (days)')
ylabel('Change from baseline of serum BCMA')
legend(labs)
set(gca,'fontsize',f.gca)
grid on

AddLetters2Plots(figure(1),{'A','B'},'HShift', -0.06, 'VShift', -0.06, ...
                'fontsize', 22)

%% Transgene copies
figure(3)
clf
hold on
TransC=0.002; % from model code
ylabel('Transgene copies/\mu g genomic DNA')
xlim([0,350])
ylim([10^0, 10^8])
set(gca,'fontsize',f.gca, 'Yscale', 'log')
grid on

% dat1
CARTe_PB = dat1.y_treat(:,1);
CARTm_PB = dat1.y_treat(:,2);
FinalCARTPB = TransC * (CARTe_PB + CARTm_PB);
plot(dat1.t_treat, FinalCARTPB, 'linewidth',lw,'color',cmap(1,:))
% dat1b
CARTe_PB = dat1b.y_treat(:,1);
CARTm_PB = dat1b.y_treat(:,2);
FinalCARTPB = TransC * (CARTe_PB + CARTm_PB);
plot(dat1b.t_treat, FinalCARTPB, 'linewidth',lw,'color',cmap(1,:),'linestyle',':')
% dat2
CARTe_PB = dat2.y_treat(:,1);
CARTm_PB = dat2.y_treat(:,2);
FinalCARTPB = TransC * (CARTe_PB + CARTm_PB);
plot(dat2.t_treat, FinalCARTPB, 'linewidth',lw,'color',cmap(2,:))
% dat3
CARTe_PB = dat3.y_treat(:,1);
CARTm_PB = dat3.y_treat(:,2);
FinalCARTPB = TransC * (CARTe_PB + CARTm_PB);
plot(dat3.t_treat, FinalCARTPB, 'linewidth',lw,'color',cmap(3,:))
%dat4
CARTe_PB = dat4.y_treat(:,1);
CARTm_PB = dat4.y_treat(:,2);
FinalCARTPB = TransC * (CARTe_PB + CARTm_PB);
plot(dat4.t_treat, FinalCARTPB, 'linewidth',lw,'color',cmap(4,:))
legend(labs)


