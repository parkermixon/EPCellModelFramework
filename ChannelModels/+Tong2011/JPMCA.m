classdef JPMCA < ChannelBase

    properties (Constant)
        num_gates = 0;
        gate_names = {""};
    end

    methods (Static)

        function Xg = Init()

            Xg = zeros(0, 0);

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

            % Evaluate current at current time
            J_PMCA_var_den = 1 + (G.K_mPMCA/Ca_i)^G.n_PMCA;
            J_PMCA_var = G.J_PMCA/J_PMCA_var_den;

            current_val = J_PMCA_var;

        end

        function [dXg, current_val] = Update(G, V, Ca_i, Xg)
            
            % Evaluate current at current time
            J_PMCA_var_den = 1 + (G.K_mPMCA/Ca_i)^G.n_PMCA;
            J_PMCA_var = G.J_PMCA/J_PMCA_var_den;
        
            current_val = J_PMCA_var;
        
            dXg = [];
        end
    end
end