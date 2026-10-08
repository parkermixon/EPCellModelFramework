classdef IClCa < ChannelBase

    properties (Constant)
        num_gates = 1;
        gate_names = {"c"};
    end

    methods (Static)

        function Xg = Init()

            Xg = zeros(1, 1);

            Xg(1) = 0.0003764413740731269; % c [unitless]

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            c = Xg(1);

            % get dXg for gating variables
            K_1Cl = 0.0006*exp((2.53*G.F*V)/(G.R*G.T));
            K_2Cl = 0.1*exp((-5*G.F*V)/(G.R*G.T));
            c_ss_den = (K_1Cl/Ca_i)^2 + (K_1Cl/Ca_i) + 1;
            c_ss = 1/(1+(K_2Cl*c_ss_den));
            tau_c_den1 = exp((V+4.56)/11.62);
            tau_c_den2 = exp((-(V+25.5))/11.62);
            tau_c = (210/(1+tau_c_den1)) + (170/(1+tau_c_den2)) - 160;
            
            dc_dt = (c_ss-c)/tau_c;
            dXg(1) = dc_dt;

        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            c = Xg(1);

            RTdivF = (G.R*G.T)/G.F;

            % Evaluate current at current time
            E_Cl = RTdivF*log(G.Cl_i/G.Cl_o);
            I_ClCa = G.g_Cl*c*(V-E_Cl);

            current_val = I_ClCa;

        end

        function [dXg, current_val] = Update(G, V, Ca_i, Xg)
            
            % Make local copies of gating variables
            c = Xg(1);

            RTdivF = (G.R*G.T)/G.F;
        
            % Evaluate current at current time
            E_Cl = RTdivF*log(G.Cl_i/G.Cl_o);
            I_ClCa = G.g_Cl*c*(V-E_Cl);
        
            current_val = I_ClCa;
        
            % get dXg for gating variables
            K_1Cl = 0.0006*exp((2.53*G.F*V)/(G.R*G.T));
            K_2Cl = 0.1*exp((-5*G.F*V)/(G.R*G.T));
            c_ss_den = (K_1Cl/Ca_i)^2 + (K_1Cl/Ca_i) + 1;
            c_ss = 1/(1+(K_2Cl*c_ss_den));
            tau_c_den1 = exp((V+4.56)/11.62);
            tau_c_den2 = exp((-(V+25.5))/11.62);
            tau_c = (210/(1+tau_c_den1)) + (170/(1+tau_c_den2)) - 160;
            dc_dt = (c_ss-c)/tau_c;
            dXg(1) = dc_dt;
        end
    end
end