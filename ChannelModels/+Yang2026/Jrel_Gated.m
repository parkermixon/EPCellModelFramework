classdef Jrel_Gated < ChannelBase

    properties (Constant)
        num_gates = 1;
        gate_names = {"Jrel_h"};
    end

    methods (Static)

        function Xg = Init()

            Xg = 0.892033577223919; % Not given, evaluated from Shuai Jang model

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]
            IP3 = X(3); % [nM]

            h = Xg(1);

            % Istim = stim(1);
            % IP3 = stim(2);

            IP3_mM = IP3 * 1e-6;

            a = 0.2; % [1/(uM*s) -> 1/(mM*ms)] equal

            alpha_h = a * G.d_2 * (IP3_mM+G.d_1)/(IP3_mM+G.d_3);
            beta_h = a*Ca_i;

            % Original
            dXg = alpha_h*(1-h) - beta_h*h;
        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]
            IP3 = X(3); % [nM]

            h = Xg(1);

            % Istim = stim(1);
            % IP3 = stim(2);

            IP3_mM = IP3 * 1e-6;

            m_ss = IP3_mM / (IP3_mM + G.d_4);
            n_ss = Ca_i / (Ca_i + G.d_5);
            
            % Eq from limit of large diffusion coefficients
            JChannel = G.v_IP3R * m_ss^3 * n_ss^3 * h^3 * (G.Ca_u - Ca_i);
            current_val = JChannel;

        end
    end

end