classdef Ib < ChannelBase

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
            
            RTdivF = (G.R*G.T)/G.F;

            % Evaluate current at current time
            E_K = RTdivF * log(G.K_o/G.K_i);
            I_b = G.g_b*(V-E_K);

            current_val = I_b;

        end

        function [dXg, current_val] = Update(G, V, Ca_i, Xg)

            RTdivF = (G.R*G.T)/G.F;
            
            % Evaluate current at current time
            E_K = RTdivF * log(G.K_o/G.K_i);
            I_b = G.g_b*(V-E_K);
        
            current_val = I_b;

            dXg = [];
        
        end
    end
end