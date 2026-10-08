classdef IKCNQ1 < ChannelBase

    properties (Constant)
        num_gates = 4;
        gate_names = {"n_Q1f", "n_Q1s", "w_Q1", "s_Q1"};
    end

    methods (Static)

        function Xg = Init()

            Xg = zeros(1, Tong2014.IKCNQ1.num_gates);

            Xg(1) = 0.090437; % n_Q1f
            Xg(2) = 0.090437; % n_Q1s
            Xg(3) = 0.923051; % w_Q1
            Xg(4) = 0.742674; % s_Q1

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [uM]

            n_Q1f = Xg(1);
            n_Q1s = Xg(2);
            w_Q1 = Xg(3);
            s_Q1 = Xg(4);

            % get dXg for gating variables
            n_Q1_ss = 1 / (1 + exp(-(V+22)/12.48));
            w_Q1_ss = 0.49 + (0.51 / (1 + exp((V+1.084)/28.78)));
            s_Q1_ss = 0.34 + (0.66 / (1 + exp((V+45.3)/12.3)));

            tau_nQ1f_den = (V + 38.1) / 33.59;
            tau_nQ1f = 395.3 / (1 + tau_nQ1f_den^2);

            tau_nQ1s_den1 = 10^(-0.02827 * (-23.9 - V));
            tau_nQ1s_den2 = 10^(-0.0357 * (V + 14.15));
            tau_nQ1s = 5503 + (5345.4 / (1 + tau_nQ1s_den1)) - (4590.6 / (1 + tau_nQ1s_den2)); % - -> + from OpenCOR

            tau_wQ1_den = (V + 48.09) / 48.83;
            tau_wQ1 = 5.44 + (29.2 / (1 + tau_wQ1_den^2));

            tau_sQ1 = 50000; % [ms]

            dnQ1f_dt = (n_Q1_ss - n_Q1f) / tau_nQ1f;
            dXg(1) = dnQ1f_dt;
            dnQ1s_dt = (n_Q1_ss - n_Q1s) / tau_nQ1s;
            dXg(2) = dnQ1s_dt;
            dwQ1_dt = (w_Q1_ss - w_Q1) / tau_wQ1;
            dXg(3) = dwQ1_dt;
            dsQ1_dt = (s_Q1_ss - s_Q1) / tau_sQ1;
            dXg(4) = dsQ1_dt;

        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [uM]

            n_Q1f = Xg(1);
            n_Q1s = Xg(2);
            w_Q1 = Xg(3);
            s_Q1 = Xg(4);

            RTdivF = (G.R*G.T)/G.F;

            % Evaluate current at current time
            E_K = RTdivF * log(G.K_o/G.K_i);
            I_KCNQ1 = G.g_KCNQ1 * (0.3 * n_Q1f + 0.7 * n_Q1s) * w_Q1 * s_Q1 * (V - E_K);

            current_val = I_KCNQ1;

        end

    end



end