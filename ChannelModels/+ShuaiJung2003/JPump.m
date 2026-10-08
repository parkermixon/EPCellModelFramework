classdef JPump < ChannelBase

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
            Ca_i = X(1); % [mM]

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
            JPump = v_P * Ca_i^2 / (k^2 + Ca_i^2);
            current_val = JPump;

        end
    end

end