classdef INaK < ChannelBase

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
            f_NaK_den = 1+0.1245*exp(-0.1*V*G.F/(G.R*G.T))+0.00219*exp(G.Na_o/49.71)*exp(-1.9*V*G.F/(G.R*G.T));
            f_NaK = 1/f_NaK_den;
            k_NaK_den = 1+((G.K_mK/G.K_o)^G.n_K);
            k_NaK = 1/k_NaK_den;
            n_NaK_den = 1+((G.K_mNa/G.Na_i)^G.n_Na);
            n_NaK = 1/n_NaK_den;
            I_NaK = G.g_NaK*f_NaK*k_NaK*n_NaK;

            current_val = I_NaK;

        end

        function [dXg, current_val] = Update(G, V, Ca_i, Xg)
            
            % Evaluate current at current time
            f_NaK_den = 1+0.1245*exp(-0.1*V*G.F/(G.R*G.T))+0.00219*exp(G.Na_o/49.71)*exp(-1.9*V*G.F/(G.R*G.T));
            f_NaK = 1/f_NaK_den;
            k_NaK_den = 1+((G.K_mK/G.K_o)^G.n_K);
            k_NaK = 1/k_NaK_den;
            n_NaK_den = 1+((G.K_mNa/G.Na_i)^G.n_Na);
            n_NaK = 1/n_NaK_den;
            I_NaK = G.g_NaK*f_NaK*k_NaK*n_NaK;
        
            current_val = I_NaK;

            dXg = [];
        
        end
    end
end