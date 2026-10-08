classdef IK2 < ChannelBase

    properties (Constant)
        num_gates = 3;
        gate_names = {"p", "k_1", "k_2"};
    end

    methods (Static)

        function Xg = Init()

            Xg = zeros(1, 3);

            Xg(1) = 0.1174074734567931; % p [unitless]
            Xg(2) = 0.9968385770271651; % k_1 [unitless]
            Xg(3) = 0.9968408069904307; % k_2 [unitless]

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            p = Xg(1);
            k_1 = Xg(2);
            k_2 = Xg(3);

            % get dXg for gating variables
            p_ss_den = 1+exp(-(V+17.91)/18.4);
            p_ss = 0.948/p_ss_den;
            k_ss_den = 1+exp((V+21.2)/5.7);
            k_ss = 1/k_ss_den;
            tau_p_den = 1+((V+64.1)/28.67)^2;
            tau_p = 100/tau_p_den;
            tau_k1_den1 = 1+exp((V-315)/50);
            tau_k1_den2 = 1+exp((-(V+74.9))/8);
            tau_k1 = 1e6*(1 - (1/(tau_k1_den1 * tau_k1_den2)));
            tau_k2_den1 = 1+exp((V-132.868)/25.3992); % CellML: 132.868, 25.3992; Appdx S1: 132.87, 25.4
            tau_k2_den2 = 1+exp((-(V+24.9203))/2.67915); % CellML: 24.9203, 2.67915; Appdx S1: 24.92, 2.68
            tau_k2 = 2.5e6*(1-(1/(tau_k2_den1*tau_k2_den2)));

            dp_dt = (p_ss-p)/tau_p;
            dXg(1) = dp_dt;
            dk1_dt = (k_ss-k_1)/tau_k1;
            dXg(2) = dk1_dt;
            dk2_dt = (k_ss-k_2)/tau_k2;
            dXg(3) = dk2_dt;

        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            p = Xg(1);
            k_1 = Xg(2);
            k_2 = Xg(3);

            RTdivF = (G.R*G.T)/G.F;

            % Evaluate current at current time
            E_K = RTdivF * log(G.K_o/G.K_i);
            I_K2 = G.g_K2*p^2*(0.75*k_1 + 0.25*k_2)*(V-E_K);

            current_val = I_K2;

        end

        function [dXg, current_val] = Update(G, V, Ca_i, Xg)
            
            % Make local copies of gating variables
            p = Xg(1);
            k_1 = Xg(2);
            k_2 = Xg(3);

            RTdivF = (G.R*G.T)/G.F;
        
            % Evaluate current at current time
            E_K = RTdivF * log(G.K_o/G.K_i);
            I_K2 = G.g_K2*p^2*(0.75*k_1 + 0.25*k_2)*(V-E_K);
        
            current_val = I_K2;
        
            % get dXg for gating variables
            p_ss_den = 1+exp(-(V+17.91)/18.4);
            p_ss = 0.948/p_ss_den;
            k_ss_den = 1+exp((V+21.2)/5.7);
            k_ss = 1/k_ss_den;
            tau_p_den = 1+((V+64.1)/28.67)^2;
            tau_p = 100/tau_p_den;
            tau_k1_den1 = 1+exp((V-315)/50);
            tau_k1_den2 = 1+exp((-(V+74.9))/8);
            tau_k1 = 1e6*(1 - (1/(tau_k1_den1 * tau_k1_den2)));
            tau_k2_den1 = 1+exp((V-132.868)/25.3992); % CellML: 132.868, 25.3992; Appdx S1: 132.87, 25.4
            tau_k2_den2 = 1+exp((-(V+24.9203))/2.67915); % CellML: 24.9203, 2.67915; Appdx S1: 24.92, 2.68
            tau_k2 = 2.5e6*(1-(1/(tau_k2_den1*tau_k2_den2)));
        
            dp_dt = (p_ss-p)/tau_p;
            dXg(1) = dp_dt;
            dk1_dt = (k_ss-k_1)/tau_k1;
            dXg(2) = dk1_dt;
            dk2_dt = (k_ss-k_2)/tau_k2;
            dXg(3) = dk2_dt;
        end
    end
end