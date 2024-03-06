function dydt = modeqns_biomarkers(t,y,Tumor, params)
    % Equations for computing biorkers
    % NOTE: requires input of the tumor values

    %% variable names
    M = y(1); % serum M-protein (biomarker)
    sBCMA = y(2); % soluble BCMA (biomarker)

    %% set parameter names
    Pm = params(1);
    gamma_m = params(2);
    Km = params(3);
    Pb = params(4);
    gamma_b = params(5);
    Kb = params(6);
    Tumor0 = params(7);
    
    % Get tumor val
    %Tumor = tumor_vals(t); % TEST THIS!

    %% Model equations
    dydt = zeros(length(y),1);
    %% M
    % d(M)/dt
    dydt(1) = Pm * (Tumor/Tumor0)^gamma_m - Km * M; % from ex code

    %% sBCMA
    % d(sBCMA)/dt    
    dydt(2) = Pb*(Tumor/Tumor0)^gamma_b - Kb * sBCMA; % from ex code
end