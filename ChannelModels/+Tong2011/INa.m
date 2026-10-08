classdef INa < ChannelBase

    properties (Constant)
        num_gates = 2;
        gate_names = {"m", "h"};
    end

    methods (Static)

        function Xg = Init()

            Xg = zeros(1, 2);

            Xg(1) = 0.1253518889572223; % m [unitless]
            Xg(2) = 0.404599170710196; % h [unitless]

        end

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            m = Xg(1);
            h = Xg(2);

            % get dXg for gating variables
            m_ss_den = 1+exp((-(V+35.9584))/9.24013); % CellML: 35.9584, 9.24013; Appdx S1: 35.96, 9.24
            m_ss = 1/m_ss_den;
            h_ss_den = 1+exp((V+57)/8);
            h_ss = 1/h_ss_den;
            tau_m_den = 1+exp((V+38)/10);
            tau_m = 0.25+(7/tau_m_den);
            tau_h_den = 1+((V+47.5)/1.5)^2;
            tau_h = 0.9+(1002.85/tau_h_den);

            dm_dt = (m_ss-m)/tau_m;
            dXg(1) = dm_dt;
            dh_dt = (h_ss-h)/tau_h;
            dXg(2) = dh_dt;

        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            m = Xg(1);
            h = Xg(2);

            RTdivF = (G.R*G.T)/G.F;

            % Evaluate current at current time
            E_Na = RTdivF*log(G.Na_o/G.Na_i);
            I_Na = G.g_Na * m^3 * h * (V-E_Na);

            current_val = I_Na;

        end

        % function [dXg, current_val] = Update(G, V, Ca_i, Xg)
        % 
        %     % Make local copies of gating variables
        %     m = Xg(1);
        %     h = Xg(2);
        % 
        %     RTdivF = (G.R*G.T)/G.F;
        % 
        %     % Evaluate current at current time
        %     E_Na = RTdivF*log(G.Na_o/G.Na_i);
        %     I_Na = G.g_Na * m^3 * h * (V-E_Na);
        % 
        %     current_val = I_Na;
        % 
        %     % get dXg for gating variables
        %     m_ss_den = 1+exp((-(V+35.9584))/9.24013); % CellML: 35.9584, 9.24013; Appdx S1: 35.96, 9.24
        %     m_ss = 1/m_ss_den;
        %     h_ss_den = 1+exp((V+57)/8);
        %     h_ss = 1/h_ss_den;
        %     tau_m_den = 1+exp((V+38)/10);
        %     tau_m = 0.25+(7/tau_m_den);
        %     tau_h_den = 1+((V+47.5)/1.5)^2;
        %     tau_h = 0.9+(1002.85/tau_h_den);
        % 
        %     dm_dt = (m_ss-m)/tau_m;
        %     dXg(1) = dm_dt;
        %     dh_dt = (h_ss-h)/tau_h;
        %     dXg(2) = dh_dt;
        % end
    end
end