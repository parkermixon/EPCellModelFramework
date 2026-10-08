classdef IKa < ChannelBase

    properties (Constant)
        num_gates = 2;
        gate_names = {"s", "x"};
    end

    methods (Static)

        function Xg = Init()

            Xg = zeros(1, 2);

            Xg(1) = 0.0307583106982354; % s [unitless]
            Xg(2) = 0.08785242843398365; % x [unitless]

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            s = Xg(1);
            x = Xg(2);

            % get dXg for gating variables
            s_ss_den = 1+exp((-(V+27.79))/7.57);
            s_ss = 1/s_ss_den;
            x_ss_den = 1+exp((V+69.5)/6);
            x_ss = 0.02+(0.98/x_ss_den);
            tau_s_den = 1+((V+20.52)/35)^2;
            tau_s = 17/tau_s_den;
            tau_x_den = 1+((V+34.18)/120)^2;
            tau_x = 7.5 + (10/tau_x_den);

            ds_dt = (s_ss-s)/tau_s;
            dXg(1) = ds_dt;
            dx_dt = (x_ss-x)/tau_x;
            dXg(2) = dx_dt;

        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            s = Xg(1);
            x = Xg(2);

            RTdivF = (G.R*G.T)/G.F;

            % Evaluate current at current time
            E_K = RTdivF * log(G.K_o/G.K_i);
            I_Ka = G.g_Ka*s*x*(V-E_K);

            current_val = I_Ka;

        end

        function [dXg, current_val] = Update(G, V, Ca_i, Xg)
            
            % Make local copies of gating variables
            s = Xg(1);
            x = Xg(2);

            RTdivF = (G.R*G.T)/G.F;
        
            % Evaluate current at current time
            E_K = RTdivF * log(G.K_o/G.K_i);
            I_Ka = G.g_Ka*s*x*(V-E_K);
        
            current_val = I_Ka;
        
            % get dXg for gating variables
            s_ss_den = 1+exp((-(V+27.79))/7.57);
            s_ss = 1/s_ss_den;
            x_ss_den = 1+exp((V+69.5)/6);
            x_ss = 0.02+(0.98/x_ss_den);
            tau_s_den = 1+((V+20.52)/35)^2;
            tau_s = 17/tau_s_den;
            tau_x_den = 1+((V+34.18)/120)^2;
            tau_x = 7.5 + (10/tau_x_den);
        
            ds_dt = (s_ss-s)/tau_s;
            dXg(1) = ds_dt;
            dx_dt = (x_ss-x)/tau_x;
            dXg(2) = dx_dt;
        end
    end
end