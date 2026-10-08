classdef SignallingParams
    properties(Constant)
        % Oxytocin Binding - Table A1
        % Maximum rate of receptor internalization
        k_int_max = 1.67e-5; % [1/ms]
        % Rate for receptor phosphorylation
        k_phos = (4e-5)*1e-6; % [mM] -> [nM]
        % Rate of degradation of OXT-OXTR complex
        k_des = 4.01e-5; % [1/ms]

        % Rate of OXT-OXTR binding
        % k_on = 4e7; % [L/(mol*min)] % units corrected from Gulliver and Pliska & Jutz
        % converted to 1/(nM*ms)
        k_on = 6.67e-7; % [1/(nM*ms)]

        % Rate of OXT-OXTR dissociation
        % k_off = 0.0224; % [1/min]; % units corrected from Gulliver and Pliska & Jutz
        % converted to [1/ms]
        k_off = 3.73e-7; % [1/ms]

        % Parameter for calculation of initial OXTR concentration
        gamma = 0.475; % Not currently used

        % IP3 Production - Original ish

        % Seem to be consistently nM, ms
        PLCbeta_tot = 133.33; % [nM]
        K_d_PLCbeta = 100; % [nM]
        k_cat = 0.0569; % [pmol/(mL*ms)] % [nM/ms] I believe
        v_5 = 0.625; % [pmol/(mL*ms)]
        k_5_IP3 = 25000 % [pmol/mL]
        v_6 = 0.625; % [pmol/(mL*ms)]
        k_6_IP3 = 25000; % [pmol/(mL*ms)]
        k_7_IP3 = 200; % [pmol/(mL*ms)]

        % IP3 Production - Yang & Jafri original vals convert to nM, ms
        % PLCbeta_tot = 133.33; % [nM]
        % K_d_PLCbeta = 100; % [nM]
        % k_cat = 0.0569; % [1/ms] from text
        % 
        % v_5 = 0.625; % [uM/s] -> [nM/ms] (same)
        % k_5_IP3 = 25e3; % [uM] -> [nM]
        % v_6 = 0.75; % [uM/s] -> [nM/ms] (same)
        % k_6_IP3 = 0.2e3; % [uM] -> [nM]
        % k_7_IP3 = 1e3; % [uM] -> [nM]
    end
end