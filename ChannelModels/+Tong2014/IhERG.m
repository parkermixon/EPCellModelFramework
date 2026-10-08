classdef IhERG < ChannelBase

    properties (Constant)
        num_gates = 3;
        gate_names = {"h_n1", "h_n2", "h_s"};
    end

    methods (Static)
        
        function Xg = Init()

            Xg = zeros(1, Tong2014.IhERG.num_gates);

            Xg(1) = 0.024990; % h_n1
            Xg(2) = 0.024990; % h_n2
            Xg(3) = 0.529214; % h_s

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [uM]

            h_n1 = Xg(1);
            h_n2 = Xg(2);
            h_s = Xg(3);

            % get dXg for gating variables
            h_n_ss = 1 / (1 + exp(-(V+16)/9.5));
            h_s_ss = 1 / (1 + exp((V+48)/24));

            tau_hn1_den1 = 1 + exp(-(V+40.8489)/13.7802);
            tau_hn1_den2 = 1 + exp((V+20.6372)/15.1113);
            tau_hn1 = 46.0999 + 1685.76 / (tau_hn1_den1*tau_hn1_den2); % Corrected from 446.09 -> 46.0999 OpenCOR

            tau_hn2_den1 = 1 + exp(-(V+41.8328)/6.96673);
            tau_hn2_den2 = 1 + exp((V+23.2432)/21.2949);
            tau_hn2 = 475.667 + 16321.6 / (tau_hn2_den1*tau_hn2_den2);

            tau_hs_den = (V + 20.7136)/44.2868;
            tau_hs = 19.7864 / (1 + tau_hs_den^2) - 0.378843;

            dhn1_dt = (h_n_ss - h_n1) / tau_hn1;
            dXg(1) = dhn1_dt;
            dhn2_dt = (h_n_ss - h_n2) / tau_hn2;
            dXg(2) = dhn2_dt;
            dhs_dt = (h_s_ss - h_s) / tau_hs;
            dXg(3) = dhs_dt;

        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [uM]

            h_n1 = Xg(1);
            h_n2 = Xg(2);
            h_s = Xg(3);

            RTdivF = (G.R*G.T)/G.F;

            % Evaluate current at current time
            E_K = RTdivF * log(G.K_o/G.K_i);
            I_hERG = G.g_hERG * (0.8 * h_n1 + 0.2 * h_n2) * h_s * (V-E_K);
            
            current_val = I_hERG;

        end

    end



end