classdef Ih

    properties (Constant)
        num_gates = 1;
        gate_names = {"y"};
    end

    methods (Static)

        function Xg = Init()

            Xg = zeros(1, 1);

            Xg(1) = 0.002604864867063448; % y [unitless]

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            y = Xg(1);

            % get dXg for gating variables
            y_ss_den = 1+exp((V+105.39)/8.6553); % CellML: 8.6553; Appdx S1: 8.66
            y_ss = 1/y_ss_den;
            tau_y_den = (3.5e-6*exp(-0.0497*V)) + (0.04003*exp(0.05211*V));  % CellML: 0.04003, 0.05211; Appdx S1: 0.04, 0.0521
            tau_y = 1/tau_y_den;

            dy_dt = (y_ss-y)/tau_y;
            dXg(1) = dy_dt;

        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            y = Xg(1);

            RTdivF = (G.R*G.T)/G.F;

            % Evaluate current at current time
            E_h_num = G.K_o + (G.P_NaK*G.Na_o);
            E_h_den = G.K_i + (G.P_NaK*G.Na_i);
            E_h = RTdivF*log(E_h_num/E_h_den);
            I_h = G.g_h * y * (V-E_h);

            current_val = I_h;

        end

        function [dXg, current_val] = Update(G, V, Ca_i, Xg)
            
            % Make local copies of gating variables
            y = Xg(1);

            RTdivF = (G.R*G.T)/G.F;
        
            % Evaluate current at current time
            E_h_num = G.K_o + (G.P_NaK*G.Na_o);
            E_h_den = G.K_i + (G.P_NaK*G.Na_i);
            E_h = RTdivF*log(E_h_num/E_h_den);
            I_h = G.g_h * y * (V-E_h);
        
            current_val = I_h;
        
            % get dXg for gating variables
            y_ss_den = 1+exp((V+105.39)/8.6553); % CellML: 8.6553; Appdx S1: 8.66
            y_ss = 1/y_ss_den;
            tau_y_den = (3.5e-6*exp(-0.0497*V)) + (0.04003*exp(0.05211*V));  % CellML: 0.04003, 0.05211; Appdx S1: 0.04, 0.0521
            tau_y = 1/tau_y_den;
        
            dy_dt = (y_ss-y)/tau_y;
            dXg(1) = dy_dt;
        end
    end
end