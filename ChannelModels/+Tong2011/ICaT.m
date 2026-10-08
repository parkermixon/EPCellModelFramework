classdef ICaT < ChannelBase

    properties (Constant)
        num_gates = 2;
        gate_names = {"b", "g"};
    end

    methods (Static)

        function Xg = Init()

            Xg = zeros(1, 2);

            Xg(1) = 0.508117603077852; % b [unitless]
            Xg(2) = 0.03582573962705717; % g [unitless]

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            b = Xg(1);
            g = Xg(2);

            % get dXg for gating variables
            b_ss_den = 1+exp((-(V+54.23))/9.88);
            b_ss = 1/b_ss_den;
            g_ss_den = 1+exp((V+72.978)/4.64);
            g_ss = 0.02 + 0.98/g_ss_den;
            tau_b_den = 1+((V+66)/26)^2;
            tau_b = 0.45 + 3.9/tau_b_den;
            tau_g_den1 = 1+exp((V-417.43)/203.18);
            tau_g_den2 = 1+exp((-(V+61.11))/8.07);
            tau_g = 150 - 150/(tau_g_den1*tau_g_den2);

            db_dt = (b_ss-b)/tau_b;
            dXg(1) = db_dt;
            dg_dt = (g_ss-g)/tau_g;
            dXg(2) = dg_dt;

        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            b = Xg(1);
            g = Xg(2);

            % Evaluate current at current time
            I_CaT = G.g_CaT * b^2 * g * (V-G.E_CaT);

            current_val = I_CaT;

        end

        function [dXg, current_val] = Update(G, V, Ca_i, Xg)
        
            % Make local copies of gating variables
            b = Xg(1);
            g = Xg(2);
        
            % Evaluate current at current time
            I_CaT = G.g_CaT * b^2 * g * (V-G.E_CaT);
        
            current_val = I_CaT;
        
            % get dXg for gating variables
            b_ss_den = 1+exp((-(V+54.23))/9.88);
            b_ss = 1/b_ss_den;
            g_ss_den = 1+exp((V+72.978)/4.64);
            g_ss = 0.02 + 0.98/g_ss_den;
            tau_b_den = 1+((V+66)/26)^2;
            tau_b = 0.45 + 3.9/tau_b_den;
            tau_g_den1 = 1+exp((V-417.43)/203.18);
            tau_g_den2 = 1+exp((-(V+61.11))/8.07);
            tau_g = 150 - 150/(tau_g_den1*tau_g_den2);
        
            db_dt = (b_ss-b)/tau_b;
            dXg(1) = db_dt;
            dg_dt = (g_ss-g)/tau_g;
            dXg(2) = dg_dt;
        end
    end
end