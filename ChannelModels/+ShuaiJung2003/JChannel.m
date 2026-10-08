classdef JChannel < ChannelBase

    properties (Constant)
        num_gates = 1;
        gate_names = {"h"};
    end

    methods (Static)

        function Xg = Init()

            Xg = 0.892033577223919; % Not given, evaluated

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            Ca_i = X(1); % [mM]

            h = Xg(1);

            IP3 = stim(1); % [concentration ?]

            % Constants
            v_C = 0.6; % [1/s]
            v_P = 0.5; % [1/s]
            v_L = 0.001; % [uM/s]
            Ca_ER = 15.0; % [uM]
            k = 0.1; % [uM]
            a = 0.2; % [1/(uM*s)]
            d_1 = 0.13; % [uM]
            d_2 = 1.05; % [uM]
            d_3 = 0.94; % [uM]
            d_4 = 0.13; % [uM]
            d_5 = 0.08; % [uM]

            % Equations
            alpha_h = a * d_2 * (IP3+d_1)/(IP3+d_3);
            beta_h = a*Ca_i;

            % Original
            dXg = alpha_h*(1-h) - beta_h*h;

            % Yang et al. 2026 mod
            % dXg = 0;
        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            Ca_i = X(1); % [mM]

            h = Xg(1);

            IP3 = stim(1); % [concentration ?]

            % Constants
            v_C = 0.6; % [1/s]
            v_P = 0.5; % [1/s]
            v_L = 0.001; % [uM/s]
            Ca_ER = 15.0; % [uM]
            k = 0.1; % [uM]
            a = 0.2; % [1/(uM*s)]
            d_1 = 0.13; % [uM]
            d_2 = 1.05; % [uM]
            d_3 = 0.94; % [uM]
            d_4 = 0.13; % [uM]
            d_5 = 0.08; % [uM]

            % Equations
            
            % Yang et al. 2026 mod
            % alpha_h = a * d_2 * (IP3+d_1)/(IP3+d_3);
            % beta_h = a*Ca_i;
            % h = alpha_h / (alpha_h + beta_h);

            m_ss = IP3 / (IP3 + d_4);
            n_ss = Ca_i / (Ca_i + d_5);
            
            % Eq from limit of large diffusion coefficients
            JChannel = v_C * m_ss^3 * n_ss^3 * h^3 * (Ca_ER - Ca_i);
            current_val = JChannel;

        end
    end

end