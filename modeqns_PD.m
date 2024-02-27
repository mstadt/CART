function dydt = modeqns_PD(t,y,params,varargin)
% PD model equations

%% variable names
Cmplx = y(1); % number of CAR-Target complexes per tumor cells (unitless?)
Nt    = y(2); % number of tumor cells (cells) in a media volume
Ne    = y(3); % number of CART cells in a media volume

%% set parameter names
MV = params(1);
K_on_CAR = params(2);
K_off_CAR = params(3);
Ag_CAR = params(4);
Ag_tumor = params(5);
Kkill_max = params(6);
KC50_CART = params(7);
DT_tumor = params(8);
DT_CART = params(9);

%% get variable inputs
% TODO

%% model equations
dydt = zeros(length(y),1);

% Cmplx
% d(Cmplx)/dt
dydt(1) = K_on_CAR * (Ag_CAR - Cmplx) * (Ag_tumor - Cmplx) - K_off_CAR * Cmplx;

% Tumor cells

SF = 1e9 / (6.023e23); % scaling factor
CplxCell = (Cmplx * MV) / (SF * Nt);
Kill = (Kkill_max * CplxCell)/(KC50_CART + CplxCell);
% d(Nt)/dt
dydt(2) = log(2)/DT_tumor * Nt - Kill * Nt;

% CAR-T cells
% d(Ne)/dt
dydt(3) = log(2)/DT_CART * Ne;

end % end modeqns