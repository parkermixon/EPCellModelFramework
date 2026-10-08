classdef IKCa < ChannelBase

    properties (Constant)
        num_gates = 2;
        gate_names = {"x_a", "x_ab1"};
    end

    methods (Static)

        function Xg = Init()

            Xg = zeros(1, 2);

            Xg(1) = 0.0003569126518797985; % x_a [unitless]
            Xg(2) = 0.002220456569762898; % x_ab1 [unitless]

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            x_a = Xg(1);
            x_ab1 = Xg(2);

            % get dXg for gating variables
            %BKa
            z_a_den1 = (((1000*Ca_i)+1538.29)/739.057)^2; % CellML: 739.057; Appdx S1: 739.06
            z_a_den2 = (((1000*Ca_i)-0.0630535)/0.161942)^2; % CellML: 0.0630535, 0.161942; Appdx S1: 0.063, 0.162
            z_a = (8.38/(1+z_a_den1))-(0.749/(1+z_a_den2)); % CellML: 8.38384, 0.749234; Appdx S1: 8.38, 0.749
            V_05a_den = (((1000*Ca_i) + 0.237503)/2.39278e-4)^0.42291; % CellML: 0.237503, ~2.39278e-4, 0.42291; Appdx S1: 0.238, 2.39e-4, 0.423
            V_05a = (5011.47/(1+V_05a_den))-37.5137; % CellML: 37.5137; Appdx S1: 37.51
            SS_a_den = exp((-z_a*G.F*(V-V_05a))/(G.R*G.T));
            SS_a = 1/(1+SS_a_den);
            tau_a_den = ((V-158.779)/(-52.1497))^2; % CellML: 158.779, 52.1497; Appdx S1: 158.78, 52.15
            tau_a = 2.40914/(1+tau_a_den); % CellML: 2.40914; Appdx S1: 2.41

            dxa_dt = (SS_a-x_a)/tau_a;
            dXg(1) = dxa_dt;

            %BKab
            z_ab1_den1 = (((1000*Ca_i)+228.71)/684.946)^2; % CellML: 684.946; Appdx S1: 684.95
            z_ab1_den2 = (((1000*Ca_i)-0.219)/0.428)^2;  % CellML: 0.218988, 0.428335
            z_ab1 = (1.4/(1+z_ab1_den1))-(0.681249/(1+z_ab1_den2)); % CellML: 0.681249; Appdx S1: 0.681
            V_05ab1_den = (((1000*Ca_i)+0.401189)/3.99115e-3)^0.668054; % CelLML: 0.401189, ~3.99115e-3, 0.668054; Appdx S1: 0.401, 3.99e-3, 0.668
            V_05ab1 = (8540.23/(1+V_05ab1_den))-109.275; % CellML: 109.275; Appdx S1: 109.28
            SS_ab1_den = exp((-z_ab1*G.F*(V-V_05ab1))/(G.R*G.T));
            SS_ab1 = 1/(1+SS_ab1_den);
            tau_ab1_den = ((V-153.019)/66.4952)^2; % CellML: 153.019, 66.4952; Appdx S1: 153.02, 66.5
            tau_ab1 = 13.8049/(1+tau_ab1_den); % CellML: 13.8049; Appdx S1: 13.8

            dxab1_dt = (SS_ab1-x_ab1)/tau_ab1;
            dXg(2) = dxab1_dt;

        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            x_a = Xg(1);
            x_ab1 = Xg(2);

            RTdivF = (G.R*G.T)/G.F;

            % Evaluate current at current time
            E_K = RTdivF * log(G.K_o/G.K_i);
            I_a = x_a * (V-E_K);
            I_ab1 = x_ab1 * (V-E_K);
            I_KCa = G.g_KCa * (G.p_a*I_a + G.p_b*I_ab1);

            current_val = I_KCa;

        end

        function [dXg, current_val] = Update(G, V, Ca_i, Xg)
            
            % Make local copies of gating variables
            x_a = Xg(1);
            x_ab1 = Xg(2);

            RTdivF = (G.R*G.T)/G.F;
        
            % Evaluate current at current time
            E_K = RTdivF * log(G.K_o/G.K_i);
            I_a = x_a * (V-E_K);
            I_ab1 = x_ab1 * (V-E_K);
            I_KCa = G.g_KCa * (G.p_a*I_a + G.p_b*I_ab1);
            
            current_val = I_KCa;
        
            % get dXg for gating variables
            %BKa
            z_a_den1 = (((1000*Ca_i)+1538.29)/739.057)^2; % CellML: 739.057; Appdx S1: 739.06
            z_a_den2 = (((1000*Ca_i)-0.0630535)/0.161942)^2; % CellML: 0.0630535, 0.161942; Appdx S1: 0.063, 0.162
            z_a = (8.38/(1+z_a_den1))-(0.749/(1+z_a_den2)); % CellML: 8.38384, 0.749234; Appdx S1: 8.38, 0.749
            V_05a_den = (((1000*Ca_i) + 0.237503)/2.39278e-4)^0.42291; % CellML: 0.237503, ~2.39278e-4, 0.42291; Appdx S1: 0.238, 2.39e-4, 0.423
            V_05a = (5011.47/(1+V_05a_den))-37.5137; % CellML: 37.5137; Appdx S1: 37.51
            SS_a_den = exp((-z_a*G.F*(V-V_05a))/(G.R*G.T));
            SS_a = 1/(1+SS_a_den);
            tau_a_den = ((V-158.779)/(-52.1497))^2; % CellML: 158.779, 52.1497; Appdx S1: 158.78, 52.15
            tau_a = 2.40914/(1+tau_a_den); % CellML: 2.40914; Appdx S1: 2.41
            
            dxa_dt = (SS_a-x_a)/tau_a;
            dXg(1) = dxa_dt;
            
            %BKab
            z_ab1_den1 = (((1000*Ca_i)+228.71)/684.946)^2; % CellML: 684.946; Appdx S1: 684.95
            z_ab1_den2 = (((1000*Ca_i)-0.219)/0.428)^2;  % CellML: 0.218988, 0.428335
            z_ab1 = (1.4/(1+z_ab1_den1))-(0.681249/(1+z_ab1_den2)); % CellML: 0.681249; Appdx S1: 0.681
            V_05ab1_den = (((1000*Ca_i)+0.401189)/3.99115e-3)^0.668054; % CelLML: 0.401189, ~3.99115e-3, 0.668054; Appdx S1: 0.401, 3.99e-3, 0.668
            V_05ab1 = (8540.23/(1+V_05ab1_den))-109.275; % CellML: 109.275; Appdx S1: 109.28
            SS_ab1_den = exp((-z_ab1*G.F*(V-V_05ab1))/(G.R*G.T));
            SS_ab1 = 1/(1+SS_ab1_den);
            tau_ab1_den = ((V-153.019)/66.4952)^2; % CellML: 153.019, 66.4952; Appdx S1: 153.02, 66.5
            tau_ab1 = 13.8049/(1+tau_ab1_den); % CellML: 13.8049; Appdx S1: 13.8
        
            dxab1_dt = (SS_ab1-x_ab1)/tau_ab1;
            dXg(2) = dxab1_dt;
        end
    end
end