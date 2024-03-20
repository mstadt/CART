% Run this after a sim to get the algebraic equations
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

CARTe_PB = y(:,1); % CARTe in blood
CARTm_PB = y(:,2); % CARTm in blood
CARTe_T  = y(:,3); % CAR-Target complexes
CARTm_T  = y(:,4); % CARTm in tumor
Cplx_T  = y(:,5); % CARTm in tumor
Tumor    = y(:,6); % tumor size

CARTe_T2PB =K21.*Vt.*CARTe_T;
CARTe_PB2T = K12.*Vb.*CARTe_PB;
CARTe_PB_deg = Kel_e.*CARTe_PB;

CARTm_T2PB = K21*Vt*CARTm_T;
CARTm_PB2T = K12*Vb*CARTm_PB;
CARTm_PB_deg = Kel_m*CARTm_PB;

CplxPCART_T = Cplx_T ./ (CARTe_T + CARTm_T);
frac_Kexp = CplxPCART_T ./ (EC50 + CplxPCART_T);
Kexp = (Kexp_max .* CplxPCART_T) ./ (EC50 + CplxPCART_T);

CARTe_exp = Kexp .* CARTe_T; % CART cell expansion
CARTe2CARTm = Rm .* CARTe_T; % conversion from effector to memory T-cells

f_CART = (CARTe_T + CARTm_T) * Ag_CAR - Cplx_T;
f_Ag   = Tumor * Ag_TAA - Cplx_T;
Cplx_Inc = Kon .* f_CART .* f_Ag;
Cplx_Dec = Koff .* Cplx_T;


CplxTumor = Cplx_T./Tumor;
frac_Kkill = CplxTumor ./ (IC50 + CplxTumor);
Kkill = (Kkill_max .* CplxTumor) ./ (IC50 + CplxTumor);

CellsKill = Kkill .* Tumor;