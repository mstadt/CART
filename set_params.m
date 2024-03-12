function p = set_params(mod_type)
if strcmp(mod_type, 'PD')
    % PD model parameters
    p.MV = 1e-6; % media volume in L (GUESSED)
    p.K_on_CAR = 7.1e4*(60*60*24)/(6.023e23); % (1/Molar/sec conv to 1/#/L/day) binding affinity of CAR to TAA (Tab 1, PD mod)
    p.K_off_CAR = 2.39e-3*(60*60*24); % (1/sec convert to 1/day) dissociation rate of CAR to TAA (Tab 1, PD mod)
    p.Ag_CAR = 15000; % numbers/CARTcell, overall density of CARs on CAR-T cells (Tab 1, PD mod)
    p.Ag_tumor = 222; % JeKo value (TODO: try others), (Tab 1, PD mod)
    p.Kkill_max = 0.353; % 1/hour # maximum rate of killing of tumor cells by CAR-T cells (Tab 1, PD mod)
    p.KC50_CART = 2.24; % (number/cell) (Tab 1, PD mod)
    p.DT_tumor = 26; % (hour), doubling time of tumor cells (Tab 1, PD mod, JeKo value)
    p.DT_CART = 24; % (hour), double time of CART cells
elseif strcmp(mod_type, 'PKPD_preclin')
    p.K12         = 20304; % (1/day) distribution rate from blood to bone marrow compartment
    p.Vb          = 0.944; % (mL) blood volume
    p.K21         = 0.3288; % (1/day) redistribution rate from bone marrow to blood compartment
    p.Vt          = 0.151; % (mL) volume of tumor compartment (or bone marrow?)
    p.Kel_e       = 113; % (1/day) elimination rate of effector CARTe (tab 1)
    p.Kel_m       = 0.219; % not given in pre-clin PKPD....this is from other part...
    p.Kexp_max    = 0.9168; % (1/day) max rate of CART cells expansion
    p.EC50_exp    = 1.15; % (num/day) 50 max for maximum rate of CART cell expansion
    p.Rm          = 0.00002; % guess from Tab 1...
    p.Ag_CAR      = 15000; % Table 1
    p.Ag_TAA      = 12590; % Table 1
    p.Kkill_max   = 0.612; % Table 1
    p.KC50_Kill   = 2.24; % Table 1
    p.Kon_CAR     = 7.1e4 * (60*60*24)/(6.023e23); % table 1 (convert to 1/#/L/day)
    p.Koff_CAR    = 2.39e-2*(60*60*24); % table 1 (convert to 1/day)
    p.Kg_tumor    = 0.0888; % first order growth rate
elseif strcmp(mod_type, 'PKPD_clin')
    p.K12         = 1.71; %20304, # (1/day) distribution rate from blood to bone marrow compartment
    p.Vb          = 5; % (L) #0.944, # (mL) blood volume
    p.K21         = 0.176; %0.3288, # (1/day) redistribution rate from bone marrow to blood compartment
    p.Vt          = 3.65; % (L) #0.151, # (mL) volume of tumor compartment (or bone marrow?)
    p.Kel_e       = 113; % (1/day) elimination rate of effector CARTe (tab 1)
    p.Kel_m       = 0.219; % not given in pre-clin PKPD....this is from other part...
    p.Kexp_max    = 1.73; %0.9168, # (1/day) max rate of CART cells expansion
    p.EC50_exp    = 10; %1.15, # (num/day) 50 max for maximum rate of CART cell expansion
    p.Rm          = 0.00002; % # guess from Tab 1...
    p.Ag_CAR      = 15000; % Table 1
    p.Ag_TAA      = 12590; % Table 1
    p.Kkill_max   = 0.343; % 0.612, # Table 1
    p.KC50_Kill   = 2.24; % Table 1
    p.Kon_CAR     = 7.1e4 * (60*60*24)/(6.023e23); % table 1 (convert to 1/#/L/day)
    p.Koff_CAR    = 2.39e-2*(60*60*24);  % table 1 (convert to 1/day)
    p.Kg_tumor    = 0.008;  %0.0888 # first order growth rate
elseif strcmp(mod_type,'biomarkers')
    M0 = 12.1*2.5E9/0.117;
    sBCMA0 = 0.5;
    p.Pm          = 12.1*2.5E9/M0; % g/cell/day * 1e-12 (from code)
    p.gamma_m     = 0.215; % Tab 1
    p.Km          = 0.117; %0.117/M0; % from code (ddt_M)
    p.Pb          = 0.7*0.5/sBCMA0; % k*B0
    p.gamma_b     = 1; % Tab 1
    p.Kb          = 0.7/sBCMA0; % from code (note -- commented out)
else
    fprintf('what is this mod_type: %s ?\n', mod_type)
    error('mod_type not supported')
end % if strcmp
end