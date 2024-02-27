function p = set_params(mod_type)
if strcmp(mod_type, 'PD')
    p.MV = 1e-6; % media volume in L (GUESSED)
    p.K_on_CAR = 7.1e4*(60*60*24)/(6.023e23); % (1/Molar/sec conv to 1/#/L/day) binding affinity of CAR to TAA (Tab 1, PD mod)
    p.K_off_CAR = 2.39e-3*(60*60*24); % (1/sec convert to 1/day) dissociation rate of CAR to TAA (Tab 1, PD mod)
    p.Ag_CAR = 15000; % numbers/CARTcell, overall density of CARs on CAR-T cells (Tab 1, PD mod)
    p.Ag_tumor = 222; % JeKo value (TODO: try others), (Tab 1, PD mod)
    p.Kkill_max = 0.353; % 1/hour # maximum rate of killing of tumor cells by CAR-T cells (Tab 1, PD mod)
    p.KC50_CART = 2.24; % (number/cell) (Tab 1, PD mod)
    p.DT_tumor = 26; % (hour), doubling time of tumor cells (Tab 1, PD mod, JeKo value)
    p.DT_CART = 24; % (hour), double time of CART cells
else
    fprintf('what is this mod_type: %s ?\n', mod_type)
    error('mod_type not supported')
end % if strcmp
end