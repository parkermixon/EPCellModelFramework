classdef IK1 < ChannelBase

    properties (Constant)
        num_gates = 3;
        gate_names = {"q", "r_1", "r_2"};
    end

    methods (Static)

        function Xg = Init()

            Xg = zeros(1, 3);

            Xg(1) = 0.2060363247740295; % q [unitless]
            Xg(2) = 0.1922244113609531; % r_1 [unitless]
            Xg(3) = 0.1932803618375963; % r_2 [unitless]

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            q = Xg(1);
            r_1 = Xg(2);
            r_2 = Xg(3);

            % get dXg for gating variables
            q_ss_den = 1+exp((-(V+18.6736))/26.6606); % CellML: 18.6736, 26.6606; Appdx S1: 18.67, 26.66
            q_ss = 0.978613/q_ss_den; % CellML: 0.978613, Appdx S1: 1
            r_ss_den = 1+exp((V+63)/6.3);
            r_ss = 1/r_ss_den;
            tau_q_den = 1+((V+60.71)/15.79)^2;
            tau_q = 500/tau_q_den;
            tau_r1_den = 1+((V+62.7133)/35.8611)^2; % CellML: 62.7133, 35.8611; Appdx S1: 62.71, 35.86
            tau_r1 = 5e3/tau_r1_den; % CellML/src: 5e3, Appdx S1: 5e4
            tau_r2_den = 1+exp((V+22)/4);
            tau_r2 = 3e4 + (2.2e5/tau_r2_den);

            dq_dt = (q_ss-q)/tau_q;
            dXg(1) = dq_dt;
            dr1_dt = (r_ss-r_1)/tau_r1;
            dXg(2) = dr1_dt;
            dr2_dt = (r_ss-r_2)/tau_r2;
            dXg(3) = dr2_dt;

        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            q = Xg(1);
            r_1 = Xg(2);
            r_2 = Xg(3);

            RTdivF = (G.R*G.T)/G.F;

            % Evaluate current at current time
            E_K = RTdivF * log(G.K_o/G.K_i);
            I_K1 = G.g_K1*q^2*(0.38*r_1 + 0.63*r_2)*(V-E_K); % Appdx S1: 0.38, 0.62; CellML, src: 0.38, 0.63

            current_val = I_K1;

        end

        function [dXg, current_val] = Update(G, V, Ca_i, Xg)
            
            % Make local copies of gating variables
            q = Xg(1);
            r_1 = Xg(2);
            r_2 = Xg(3);

            RTdivF = (G.R*G.T)/G.F;
        
            % Evaluate current at current time
            E_K = RTdivF * log(G.K_o/G.K_i);
            I_K1 = G.g_K1*q^2*(0.38*r_1 + 0.63*r_2)*(V-E_K); % Appdx S1: 0.38, 0.62; CellML, src: 0.38, 0.63
        
            current_val = I_K1;
        
            % get dXg for gating variables
            q_ss_den = 1+exp((-(V+18.6736))/26.6606); % CellML: 18.6736, 26.6606; Appdx S1: 18.67, 26.66
            q_ss = 0.978613/q_ss_den; % CellML: 0.978613, Appdx S1: 1
            r_ss_den = 1+exp((V+63)/6.3);
            r_ss = 1/r_ss_den;
            tau_q_den = 1+((V+60.71)/15.79)^2;
            tau_q = 500/tau_q_den;
            tau_r1_den = 1+((V+62.7133)/35.8611)^2; % CellML: 62.7133, 35.8611; Appdx S1: 62.71, 35.86
            tau_r1 = 5e3/tau_r1_den; % CellML/src: 5e3, Appdx S1: 5e4
            tau_r2_den = 1+exp((V+22)/4);
            tau_r2 = 3e4 + (2.2e5/tau_r2_den);
        
            dq_dt = (q_ss-q)/tau_q;
            dXg(1) = dq_dt;
            dr1_dt = (r_ss-r_1)/tau_r1;
            dXg(2) = dr1_dt;
            dr2_dt = (r_ss-r_2)/tau_r2;
            dXg(3) = dr2_dt;
        end
    end
end