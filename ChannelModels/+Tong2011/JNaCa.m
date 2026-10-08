classdef JNaCa < ChannelBase

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
            f_allo_den = 1+((G.K_mAllo/Ca_i)^G.n_Allo);
            f_allo = 1/f_allo_den;
            f_1NaCa = exp((G.gam-1)*((V*G.F)/(G.R*G.T)));
            f_2NaCa = exp(G.gam*((V*G.F)/(G.R*G.T)));
            J_NaCa_var_den1 = 1+G.k_sat*f_1NaCa; % 1st coefficient
            J_NaCa_var_den2 = (G.K_mCao*G.Na_i^3) + ((G.K_mNao^3)*Ca_i) + (G.Ca_o*G.Na_i^3) + (G.Na_o^3*Ca_i); % 2nd coefficient; 1st half of sum
            J_NaCa_var_den3 = ((G.K_mNai^3)*G.Ca_o)*(1+(Ca_i/G.K_mCai)) + (G.Na_o^3*G.K_mCai)*(1+(G.Na_i/G.K_mNai)^3); % 2nd coefficient; 2nd half of sum
            J_NaCa_var_num1 = G.J_NaCa*f_allo;
            J_NaCa_var_num2 = (G.Na_i^3*G.Ca_o*f_2NaCa)-(G.Na_o^3*Ca_i*f_1NaCa);
            % Equations for I_NaCa/J_NaCa vary between all 3 sources; both CellML & src
            % include the 0.5 factor, while Appdx S1 does not
            J_NaCa_var = -(J_NaCa_var_num1*J_NaCa_var_num2) / (J_NaCa_var_den1*(J_NaCa_var_den2+J_NaCa_var_den3));

            current_val = J_NaCa_var;

        end

        function [dXg, current_val] = Update(G, V, Ca_i, Xg)
        
            % Evaluate current at current time
            f_allo_den = 1+((G.K_mAllo/Ca_i)^G.n_Allo);
            f_allo = 1/f_allo_den;
            f_1NaCa = exp((G.gam-1)*((V*G.F)/(G.R*G.T)));
            f_2NaCa = exp(G.gam*((V*G.F)/(G.R*G.T)));
            J_NaCa_var_den1 = 1+G.k_sat*f_1NaCa; % 1st coefficient
            J_NaCa_var_den2 = (G.K_mCao*G.Na_i^3) + ((G.K_mNao^3)*Ca_i) + (G.Ca_o*G.Na_i^3) + (G.Na_o^3*Ca_i); % 2nd coefficient; 1st half of sum
            J_NaCa_var_den3 = ((G.K_mNai^3)*G.Ca_o)*(1+(Ca_i/G.K_mCai)) + (G.Na_o^3*G.K_mCai)*(1+(G.Na_i/G.K_mNai)^3); % 2nd coefficient; 2nd half of sum
            J_NaCa_var_num1 = G.J_NaCa*f_allo;
            J_NaCa_var_num2 = (G.Na_i^3*G.Ca_o*f_2NaCa)-(G.Na_o^3*Ca_i*f_1NaCa);
            % Equations for I_NaCa/J_NaCa vary between all 3 sources; both CellML & src
            % include the 0.5 factor, while Appdx S1 does not
            J_NaCa_var = -(J_NaCa_var_num1*J_NaCa_var_num2) / (J_NaCa_var_den1*(J_NaCa_var_den2+J_NaCa_var_den3));
        
            current_val = J_NaCa_var;

            dXg = [];
        
        end
    end
end