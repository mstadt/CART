function dydt = modeqns_PKPD_CARTe_sep(t,y,params,varargin)
% PKPD model equations

%% variable names
CARTe_PB = y(1); % CARTe in blood
CARTm_PB = y(2); % CARTm in blood
CARTe_T  = y(3); % CARTe in tissue from moving
CARTm_T  = y(4); % CARTm in tissue
Cplx_T     = y(5); % CAR-Target complexes
Tumor    = y(6); % tumor size
CARTe_exp = y(7);
B = y(8);

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
Pm = params(17);
gamma_m = params(18);
Km = params(19);
Pb = params(20);
gamma_b = params(21);
Kb = params(22);
Tumor0 = params(23);
M0 = params(24);
B0 = params(25);


%% get variable inputs
% doseCART = 0; % default no treatment
% for ii = 1:2:length(varargin)
%     temp = varargin{ii+1};
%     if strcmp(varargin{ii},'doseCART')
%         doseCART = temp;
%     end
% end

%% model equations
dydt = zeros(length(y),1);

%% Peripheral blood compartment
% d(CARTe_PB)/dt
CARTe_T2PB = 0; %K21*Vt*CARTe_T; % CARTe cell T to PB compartment
CARTe_PB2T = K12*Vb*CARTe_PB; % CARTe cell PB to T compartment
CARTe_PB_deg = 0; %Kel_e*CARTe_PB; % CARTe degradation in PB
dydt(1) = (CARTe_T2PB - CARTe_PB2T)./Vb - CARTe_PB_deg;

% d(CARTm_PB)/dt
CARTm_T2PB = 0; %K21*Vt*CARTm_T; % CARTm cell T to PB compartment
CARTm_PB2T = 0; %K12*Vb*CARTm_PB; % CARTm cell PB to T compartment
CARTm_PB_deg = 0; %Kel_m*CARTm_PB; % degradation of CARTm in PB
dydt(2) = (CARTm_T2PB - CARTm_PB2T)./Vb  - CARTm_PB_deg;

%% Tissue compartment (bone marrow or solid tumor)
if (CARTe_T + CARTm_T) > 0
    CplxPCART_T = Cplx_T / (CARTe_T + CARTm_T);  
    % CAR-T cell expansion
    Kexp = (Kexp_max * CplxPCART_T) ./ (EC50_exp + CplxPCART_T);
else
    Kexp = 0;
end

% d(CARTe_T)/dt
CARTe_exp = Kexp * CARTe_T; % CART cell expansion
CARTe2CARTm = 0; %Rm * CARTe_T; % conversion from effector to memory T-cells
dydt(3) = (CARTe_PB2T - CARTe_T2PB)./Vt - CARTe2CARTm;

dydt(7) = CARTe_exp;

% d(CARTm_T)/dt
dydt(4) = (CARTm_PB2T - CARTm_T2PB)./Vt + CARTe2CARTm;

% Cplx_T
f_CART = (CARTe_T + CARTm_T) * Ag_CAR - Cplx_T;
f_Ag   = Tumor * Ag_TAA - Cplx_T;
% d(Cplx)/dt
dydt(5) = Kon_CAR * f_CART * f_Ag - Koff_CAR * Cplx_T;

%% Tumor
CplxTumor = Cplx_T./Tumor;
Kkill = (Kkill_max * CplxTumor) / (KC50_Kill + CplxTumor);

% d(Tumor)/dt
dydt(6) = Kg_tumor * Tumor - Kkill * Tumor;

% dM/dt 

%dydt(7) = 0; %Pm/M0*(Tumor/Tumor0).^gamma_m - Km*M/M0;

% dBdt
dydt(8) = 0; %Pb/B0*(Tumor/Tumor0).^gamma_b - Kb*B/B0;

end % end modeqns