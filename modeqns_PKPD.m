function dydt = modeqns_PKPD(t,y,params,varargin)
% PKPD model equations

%% variable names
CARTe_PB = y(1); % CARTe in blood
CARTm_PB = y(2); % CARTm in blood
CARTe_T  = y(3); % CARTe in tissue
CARTm_T  = y(4); % CARTm in tissue
Cplx     = y(5); % CAR-Target Complexes
Tumor    = y(6); % tumor size

%% set parameter names
K12 = params(1);
Vb = params(2);
K21 = params(3);
Vt = params(4);
Kel_e = params(5);
Kel_m = params(6);
Kexp_max = params(7);
EC50_exp = params(8);
Rm = params(9);
Ag_CAR = params(10);
Ag_TAA = params(11);
Kkill_max = params(12);
KC50_Kill = params(13);
Kon_CAR = params(14);
Koff_CAR = params(15);
Kg_tumor = params(16);

%% get variable inputs
doseCART = 0; % default no treatment
for ii = 1:2:length(varargin)
    temp = varargin{ii+1};
    if strcmp(varargin{ii},'doseCART')
        doseCART = temp;
    end
end

%% model equations
dydt = zeros(length(y),1);

%% Peripheral blood compartment
% d(CARTe_PB)/dt
dydt(1) = (doseCART - K12*Vb*CARTe_PB + K21*Vt*CARTe_T)./Vb - Kel_e*CARTe_PB;

% d(CARTm_PB)/dt
dydt(2) = (-K12*Vb*CARTm_PB + K21*Vt*CARTm_T)./Vb  - Kel_m*CARTm_PB;

%% Tissue compartment (bone marrow or solid tumor)
CplxCART = Cplx ./ (CARTe_T + CARTm_T);
Kexp = (Kexp_max * CplxCART) ./ (EC50_exp + CplxCART);
% d(CARTe_T)/dt
dydt(3) = (K12*Vb*CARTe_PB - K21*Vt*CARTe_T)./Vt + Kexp*CARTe_T - Rm*CARTe_T;

% d(CARTm_T)/dt
dydt(4) = (K12*Vb*CARTm_PB + K21*Vt*CARTm_T)./Vt + Rm*CARTe_T;

% Cplx
f_CART = (CARTe_T + CARTm_T) * Ag_CAR - Cplx;
f_Ag   = Tumor * Ag_TAA - Cplx;
% d(Cplx)/dt
dydt(5) = Kon_CAR * f_CART * f_Ag - Koff_CAR * Cplx;

%% Tumor
CplxTumor = Cplx./Tumor;
Kkill = (Kkill_max * CplxTumor) / (KC50_Kill + CplxTumor);

% d(Tumor)/dt
dydt(6) = Kg_tumor * Tumor - Kkill * CARTe_T * Tumor;


end % end modeqns