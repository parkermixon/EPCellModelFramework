classdef Jrel < ChannelBase

    properties (Constant)
        num_gates = 0;
        gate_names = {""};
    end

    methods (Static)

        function Xg = Init()

            Xg = [];

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            dXg = [];

        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]
            IP3 = X(3); % [nM]
            % Istim = stim(1);
            % IP3 = stim(2);

            IP3_mM = IP3 * 1e-6;

            % Removed a as it is cancelled out in h formula
            alpha = G.d_2 * (IP3_mM + G.d_1)/(IP3_mM + G.d_3);
            beta = Ca_i;

            h = alpha / (alpha + beta);

            prod1 = IP3_mM / (G.d_4 + IP3_mM);
            prod2 = Ca_i / (G.d_5 + Ca_i);

            Jrel = G.v_IP3R * (prod1 * prod2 * h)^3 * (G.Ca_u - Ca_i);

            current_val = Jrel;

        end
    end

end