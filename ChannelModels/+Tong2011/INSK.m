classdef INSK < ChannelBase

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
            % named g_NS_var to prevent overwritting/confusing with constant g_NS
            P_CaCs_prime = G.P_CaCs/(1+exp((V*G.F)/(G.R*G.T)));
            E_NS_num = G.P_NaCs*G.Na_o + G.P_KCs*G.K_o + 4*P_CaCs_prime*G.Ca_o;
            E_NS_den = G.P_NaCs*G.Na_i + G.P_KCs*G.K_i + 4*P_CaCs_prime*Ca_i;
            E_NS = RTdivF * log(E_NS_num/E_NS_den);

            f_Mg_den = 1 + (G.Mg_o/G.K_dMg)^1.29834; % CellML: 1.29834; Appdx S1: 1.3
            f_Mg = 0.108043 + (0.903902/f_Mg_den); % CellML: 0.108043, 0.903902; Appdx S1: 0.1, 0.9

            g_sNa = 0.0123;
            g_sK = g_sNa;
            g_of_K_den = 1+(150/(G.K_o+1e-8))^2;
            g_of_K = (1/g_sK) *(0.03/g_of_K_den);
            g_NSK = G.g_NS*1.19*g_of_K; % in CellML variable gnsK=1.19 is used instead of 1.19

            I_NSK = g_NSK*f_Mg*(V-E_NS);

            current_val = I_NSK;

        end

        function [dXg, current_val] = Update(G, V, Ca_i, Xg)

            RTdivF = (G.R*G.T)/G.F;
        
            % Evaluate current at current time
            % named g_NS_var to prevent overwritting/confusing with constant g_NS
            P_CaCs_prime = G.P_CaCs/(1+exp((V*G.F)/(G.R*G.T)));
            E_NS_num = G.P_NaCs*G.Na_o + G.P_KCs*G.K_o + 4*P_CaCs_prime*G.Ca_o;
            E_NS_den = G.P_NaCs*G.Na_i + G.P_KCs*G.K_i + 4*P_CaCs_prime*Ca_i;
            E_NS = RTdivF * log(E_NS_num/E_NS_den);
        
            f_Mg_den = 1 + (G.Mg_o/G.K_dMg)^1.29834; % CellML: 1.29834; Appdx S1: 1.3
            f_Mg = 0.108043 + (0.903902/f_Mg_den); % CellML: 0.108043, 0.903902; Appdx S1: 0.1, 0.9
            
            g_sNa = 0.0123;
            g_sK = g_sNa;
            g_of_K_den = 1+(150/(G.K_o+1e-8))^2;
            g_of_K = (1/g_sK) *(0.03/g_of_K_den);
            g_NSK = G.g_NS*1.19*g_of_K; % in CellML variable gnsK=1.19 is used instead of 1.19
        
            I_NSK = g_NSK*f_Mg*(V-E_NS);
        
            current_val = I_NSK;

            dXg = [];
        
        end
    end
end