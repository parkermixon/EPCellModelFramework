classdef IKCNQ4 < ChannelBase

    properties (Constant)
        num_gates = 2;
        gate_names = {"n_Q4", "s_Q4"};
    end

    methods (Static)

        function Xg = Init()

            Xg = zeros(1, Tong2014.IKCNQ4.num_gates);

            Xg(1) = 0.108110; % n_Q4
            Xg(2) = 0.628062; % s_Q4

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [uM]

            n_Q4 = Xg(1);
            s_Q4 = Xg(2);

            % get dXg for gating variables
            n_Q4_ss = 1 / (1 + exp(-(V+15.04)/16.95));
            
            s_Q4_den1 = exp((V+86.84)/15.05);
            s_Q4_den2 = exp((V-70.13)/13.37);
            s_Q4_ss = (0.4058 / (1 + s_Q4_den1)) + (0.5942 / (1 + s_Q4_den2)); 

            tau_nQ4_den = 1 + exp((-18.01 - V)/31.04);
            tau_nQ4 = 10 + (895.9 / tau_nQ4_den);

            tau_sQ4_den = (V-39.44)/7.344;
            tau_sQ4 = 1077 + (185845/(1+tau_sQ4_den^2));

            dnQ4_dt = (n_Q4_ss - n_Q4) / tau_nQ4;
            dXg(1) = dnQ4_dt;
            dsQ4_dt = (s_Q4_ss - s_Q4) / tau_sQ4;
            dXg(2) = dsQ4_dt;

        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [uM]

            n_Q4 = Xg(1);
            s_Q4 = Xg(2);

            RTdivF = (G.R*G.T)/G.F;

            % Evaluate current at current time
            E_K = RTdivF * log(G.K_o/G.K_i);
            I_KCNQ4 = G.g_KCNQ4 * n_Q4 * s_Q4 * (V - E_K);

            current_val = I_KCNQ4;

        end

    end



end