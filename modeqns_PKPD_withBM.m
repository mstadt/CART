function dydt = modeqns_PKPD_withBM(t,y,params,varargin)
% PKPD model equations

%% variable names
CARTe_PB   = y(1); % CARTe in blood
CARTm_PB   = y(2); % CARTm in blood
CARTe_T    = y(3); % CARTe in tissue
CARTm_T    = y(4); % CARTm in tissue
Cplx_T     = y(5); % CAR-Target complexes
Tumor_T    = y(6); % tumor size
M = y(7);
B = y(8);

%% set parameter names
K12 = params(1);
Vb = params(2);
K21 = params(3);
Vt = params(4);
Kel_e = params(5);
Kel_m = params(6);
Kexp_max = params(7);
EC50 = params(8);
Rm = params(9);
Ag_CAR = params(10);
Ag_TAA = params(11);
Kkill_max = params(12);
IC50 = params(13);
Kon = params(14);
Koff = params(15);
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

%% model equations
dydt = zeros(length(y),1);

%% Peripheral blood compartment
% d(CARTe_PB)/dt
CARTe_T2PB = K21*Vt*CARTe_T; % CARTe cell T to PB compartment
CARTe_PB2T = K12*Vb*CARTe_PB; % CARTe cell PB to T compartment
CARTe_PB_deg = Kel_e*CARTe_PB; % CARTe degradation in PB
dydt(1) = (CARTe_T2PB - CARTe_PB2T)./Vb - CARTe_PB_deg;

% d(CARTm_PB)/dt
CARTm_T2PB = K21*Vt*CARTm_T; % CARTm cell T to PB compartment
CARTm_PB2T = K12*Vb*CARTm_PB; % CARTm cell PB to T compartment
CARTm_PB_deg = Kel_m*CARTm_PB; % degradation of CARTm in PB
dydt(2) = (CARTm_T2PB - CARTm_PB2T)./Vb - CARTm_PB_deg;

%% Tissue compartment (bone marrow or solid tumor)
Density_CAR = Ag_CAR;
Density_TAA = Ag_TAA;
CAR_T = (CARTe_T+ CARTm_T) * Density_CAR - Cplx_T;  
TAA_T = Tumor_T * Density_TAA - Cplx_T;
CplxPTumor_T= Cplx_T / Tumor_T;
if (CARTe_T + CARTm_T) > 0
    CplxPCART_T= Cplx_T / (CARTe_T + CARTm_T);
else
    CplxPCART_T= 0; %Cplx_T / (1e-6 + CARTm_T);
end 

% CARTe_T
Kexp = (Kexp_max* CplxPCART_T/ (EC50 + CplxPCART_T));
dydt(3) = (CARTe_PB2T - CARTe_T2PB)/Vt + Kexp * CARTe_T - Rm*CARTe_T;

% d(CARTm_T)/dt
dydt(4) = (CARTm_PB2T - CARTm_T2PB)./Vt + Rm*CARTe_T;

% d(Cplx_T)/dt
dydt(5) = Kon * CAR_T * TAA_T - Koff * Cplx_T; 

% d(Tumor)/dt
dydt(6) = Kg_tumor * Tumor_T - (Kkill_max * (CplxPTumor_T) / (IC50 + (CplxPTumor_T))) * Tumor_T;

% dM/dt
dydt(7) = Pm*(Tumor_T/Tumor0).^gamma_m - Km*M;
% 
% dBdt
dydt(8) = Pb*(Tumor_T/Tumor0).^gamma_b - Kb*B;

end % end modeqns