classdef IKCNQ5 < ChannelBase

    properties (Constant)
        num_gates = 4;
        gate_names = {"n_Q5f", "n_Q5s", "w_Q5", "s_Q5"};
    end

    methods (Static)

        function Xg = Init()

            Xg = zeros(1, Tong2014.IKCNQ5.num_gates);

            Xg(1) = 0.261889; % n_Q5f
            Xg(2) = 0.261889; % n_Q5s
            Xg(3) = 0.923051; % w_Q5
            Xg(4) = 0.628062; % s_Q5

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [uM]

            n_Q5f = Xg(1);
            n_Q5s = Xg(2);
            w_Q5 = Xg(3);
            s_Q5 = Xg(4);

            % get dXg for gating variables

            % Steady States
            n_Q5_ss = 1 / (1 + exp(-(V+36.55)/13.76));

            % w_Q5_ss = w_Q1_ss
            w_Q5_ss = 0.49 + (0.51 / (1 + exp((V+1.084)/28.78)));

            % s_Q5_ss = s_Q4_ss
            s_Q5_den1 = exp((V+86.84)/15.05);
            s_Q5_den2 = exp((V-70.13)/13.37);
            s_Q5_ss = (0.4058 / (1 + s_Q5_den1)) + (0.5942 / (1 + s_Q5_den2)); 

            % Time Constants
            tau_nQ5f_den = (V + 40.24) / 17.72;
            tau_nQ5f = 37.51 + (539 / (1 + tau_nQ5f_den^2));

            tau_nQ5s = 1000; % [ms]

            % tau_wQ5 = tau_wQ1
            tau_wQ5_den = (V + 48.09) / 48.83;
            tau_wQ5 = 5.44 + (29.2 / (1 + tau_wQ5_den^2));

            % tau_sQ5 = tau_sQ4
            tau_sQ5_den = (V-39.44)/7.344;
            tau_sQ5 = 1077 + (185845/(1+tau_sQ5_den^2));

            dnQ5f_dt = (n_Q5_ss - n_Q5f) / tau_nQ5f;
            dXg(1) = dnQ5f_dt;
            dnQ5s_dt = (n_Q5_ss - n_Q5s) / tau_nQ5s;
            dXg(2) = dnQ5s_dt;
            dwQ5_dt = (w_Q5_ss - w_Q5) / tau_wQ5;
            dXg(3) = dwQ5_dt;
            dsQ5_dt = (s_Q5_ss - s_Q5) / tau_sQ5;
            dXg(4) = dsQ5_dt;

        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [uM]

            n_Q5f = Xg(1);
            n_Q5s = Xg(2);
            w_Q5 = Xg(3);
            s_Q5 = Xg(4);

            RTdivF = (G.R*G.T)/G.F;

            % Evaluate current at current time
            E_K = RTdivF * log(G.K_o/G.K_i);
            I_KCNQ5 = G.g_KCNQ5 * (0.2 * n_Q5f + 0.8 * n_Q5s) * w_Q5 * s_Q5 * (V - E_K);

            current_val = I_KCNQ5;

        end

    end



end