classdef OXT_OXTR < ChannelBase % TODO: implement as a sep signal class

% NOT USED

    properties (Constant)
        num_gates = 4;
        gate_names = {"C_oxt", "C_oxtr", "C_oxtrc", "C_oxtrc_int"};
    end

    methods (Static)

        function Xg = Init()

            Xg = zeros(1, 2);

            Xg(1) = 0; % C_oxt [unitless] % Should be stim
            Xg(2) = 0; % C_oxtr [unitless]
            Xg(3) = 0; % C_oxtrc
            Xg(4) = 0; % C_oxtrc_int

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]
            IP3 = X(3); % [nM]

            C_oxt = Xg(1); % TODO: Move to stim
            C_oxtr = Xg(2);
            C_oxtrc = Xg(3);
            C_oxtrc_int = Xg(4);

            % Initialize outputs to zero
            dXg = zeros(size(Xg));

            % Evaluate rate constants
            k_int = (G.k_int_max * C_oxtrc) / (G.k_phos + C_oxtrc);
            k_rec = k_int;

            % OXTR Pathways
            OXTR_OXTRC_on = G.k_on * C_oxt * C_oxtr;
            OXTR_OXTRC_off = G.k_off * C_oxtrc;
            OXTR_OXTRC_INT_rec = k_rec * C_oxtrc_int;
            OXTR_OXTRC_INT_int = k_int * C_oxtrc;

            dCoxtr_dt = -OXTR_OXTRC_on + OXTR_OXTRC_off + OXTR_OXTRC_INT_rec;
            dXg(2) = dCoxtr_dt;

            dCoxtrc_dt = OXTR_OXTRC_on - OXTR_OXTRC_off - OXTR_OXTRC_INT_int - (G.k_des * C_oxtrc);
            dXg(3) = dCoxtrc_dt;

            dCoxtrc_int_dt = OXTR_OXTRC_INT_int - OXTR_OXTRC_INT_rec;
            dXg(4) = dCoxtrc_int_dt;
        end

        % Update current value
        % Essentially just a wrapper to return C_oxtrc to model
        % Should be fixed to actually be it's own signalling stuff
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]
            IP3 = X(3); % [nM]

            C_oxt = Xg(1); % TODO: Move to stim
            % C_oxtr = Xg(2);
            C_oxtrc = Xg(3);
            % C_oxtrc_int = Xg(4);

            current_val = C_oxtrc;

        end
    end

end