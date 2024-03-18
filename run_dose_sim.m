function [t,y] = run_dose_sim(doseCART_tot, dose_start,dose_time_hrs, params, tspan, IC, options)
    % Runs a simulation that includes a dose of CART cells
    % Inputs:
    %   doseCART_tot -- total CART cells in dose
    %   dose_start -- time to start dose
    %   dose_time_hrs -- how long dose is given (in hours)
    %   params -- parameters for model
    %   tspan -- start to end time of simulation
    %   IC -- initial condition
    %   options -- ODE solver settings
    doseCART = doseCART_tot/(dose_time_hrs/24);
    if dose_start > tspan(1)
        % start with no dose to start
        t0 = tspan(1);
        tf = dose_start;
        IC1 = IC;
        [t1,y1] = ode15s(@(t,y) modeqns_PKPD(t,y,params,...
                                            'doseCART', 0),...
                                            [t0,tf], IC1, options);
        % add dose
        IC2 = y1(end,:);
        t0 = dose_start;
        tf = min(dose_start + dose_time_hrs/24, tspan(2));
        [t2,y2] = ode15s(@(t,y) modeqns_PKPD(t,y,params,...
                                            'doseCART', doseCART),...
                                            [t0,tf], IC2, options);

        t = [t1;t2];
        y = [y1;y2];
        if tspan(2) > (dose_start + dose_time_hrs/24)
            IC3 = y2(end,:);
            t0 = t2(end);
            tf = tspan(2);
            [t3,y3] = ode15s(@(t,y) modeqns_PKPD(t,y,params,...
                                            'doseCART', 0),...
                                            [t0,tf], IC3, options);
            t = [t;t3];
            y = [y;y3];
        end
    elseif dose_start == tspan(1)
        % start with dose
        IC1 = IC;
        t0 = tspan(1);
        tf = dose_time_hrs/24;
        [t1,y1] = ode15s(@(t,y) modeqns_PKPD(t,y,params,...
                                            'doseCART', doseCART),...
                                            [t0,tf], IC1, options);
        t = t1;
        y = y1;
       if tspan(2) > (dose_start + dose_time_hrs/24)
            IC2 = y1(end,:);
            t0 = t1(end);
            tf = tspan(2);
            [t2,y2] = ode15s(@(t,y) modeqns_PKPD(t,y,params,...
                                            'doseCART', 0),...
                                            [t0,tf], IC2, options);
            t = [t;t2];
            y = [y;y2];
        end
    else
        fprintf('dose_start: %f, tspan(1): %f \n', dose_start, tspan(1))
        error('dose_start before simulation start')
    end
end