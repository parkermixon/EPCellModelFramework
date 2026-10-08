classdef ICaL < ChannelBase

    properties (Constant)
        num_gates = 3;
        gate_names = {"d", "f_1", "f_2"};
    end

    methods (Static)

        function Xg = Init()
            
            Xg = zeros(1, 3);

            Xg(1) = 0.01036961357784695; % d [unitless]
            Xg(2) = 0.9065941499695301; % f_1 [unitless]
            Xg(3) = 0.9065967263076083; % f_2 [unitless]

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]
            
            d = Xg(1);
            f_1 = Xg(2);
            f_2 = Xg(3);

            % get dXg for gating variables
            d_ss_den = 1 + exp(-(V+22)/7);
            d_ss = 1/d_ss_den;
            f_ss_den = 1 + exp((V+38)/G.f_ss_sf); % inactivation steady state
            f_ss = 1/f_ss_den;
            tau_d_den = 1+((V+29.97)/9)^2;
            tau_d = 2.29 + 5.7/tau_d_den;
            tau_f1 = 12; % [ms]
            tau_f2_den1 = (1+exp((V+13.9629)/45.3782));
            tau_f2_den2 = (1+exp(-(V+9.49866)/3.3945));
            tau_f2 = 90.9699 * (1-(1/(tau_f2_den1*tau_f2_den2)));

            dd_dt = (d_ss-d)/tau_d;
            dXg(1) = dd_dt;
            df1_dt = (f_ss-f_1)/tau_f1;
            dXg(2) = df1_dt;
            df2_dt = (f_ss-f_2)/tau_f2;
            dXg(3) = df2_dt;

        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]
            
            d = Xg(1);
            f_1 = Xg(2);
            f_2 = Xg(3);

            % Evaluate current at current time
            f_Ca_den = 1 + (Ca_i/G.K_mCaL)^4;
            f_Ca = 1 / f_Ca_den;
            I_CaL = G.g_CaL * d^2 * f_Ca * (0.8*f_1 + 0.2*f_2)*(V-G.E_CaL);

            current_val = I_CaL;

        end

        function [dXg, current_val] = Update(G, V, Ca_i, Xg)
    
            % Make local copies of gating variables
            d = Xg(1);
            f_1 = Xg(2);
            f_2 = Xg(3);
            
            % Evaluate current at current time
            f_Ca_den = 1 + (Ca_i/G.K_mCaL)^4;
            f_Ca = 1 / f_Ca_den;
            I_CaL = G.g_CaL * d^2 * f_Ca * (0.8*f_1 + 0.2*f_2)*(V-G.E_CaL);
        
            current_val = I_CaL;
        
            % get dXg for gating variables
            d_ss_den = 1 + exp(-(V+22)/7);
            d_ss = 1/d_ss_den;
            f_ss_den = 1 + exp((V+38)/G.f_ss_sf); % inactivation steady state
            f_ss = 1/f_ss_den;
            tau_d_den = 1+((V+29.97)/9)^2;
            tau_d = 2.29 + 5.7/tau_d_den;
            tau_f1 = 12; % [ms]
            tau_f2_den1 = (1+exp((V+13.9629)/45.3782));
            tau_f2_den2 = (1+exp(-(V+9.49866)/3.3945));
            tau_f2 = 90.9699 * (1-(1/(tau_f2_den1*tau_f2_den2)));
            
            dd_dt = (d_ss-d)/tau_d;
            dXg(1) = dd_dt;
            df1_dt = (f_ss-f_1)/tau_f1;
            dXg(2) = df1_dt;
            df2_dt = (f_ss-f_2)/tau_f2;
            dXg(3) = df2_dt;
    
        end
    end
end