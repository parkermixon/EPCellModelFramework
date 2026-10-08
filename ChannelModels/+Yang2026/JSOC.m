classdef JSOC < ChannelBase

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
            % IP3 = stim(2); % [concentration ?]

            JSOC = G.k_SRleak * (G.Ca_u - Ca_i);

            current_val = JSOC;

        end
    end

end