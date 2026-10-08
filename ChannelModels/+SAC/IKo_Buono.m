classdef IKo_Buono %< ChannelBase
    
    % TREK-1

    % Does not account for K_o or K_i as original model does

    properties (Constant)
        num_gates = 0;
        gate_names = {""};

        % Values from Buonocunto source code (protocols.m, model_SAC.m)
        % g_Ko = 12.8208e-6;
        % beta_Ko = 0.8742;

    end

    methods (Static)
        
        function Xg = Init()

            Xg = [];
            
        end

        function dXg = GetG(G, X, Xg, stim)

            dXg = [];

        end

        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM?]

            % Istim = stim(1);
            lambda = stim(2); % stretch [%]

            % R=8314.0;
            % T=310.0;
            % F=96485.0;
            % 
            % ko=5.4; %only for I-V curves
            % ki=142.412524737626; % Initial cond, using as constant

            E_K=(G.R*G.T/G.F)*log(G.ko/G.ki);

            eps = lambda - 1;

            gamma_SL_Ko = G.beta_Ko * eps + 0.7;

            I_Ko = G.g_Ko * (gamma_SL_Ko / (1 + exp(-((10 + V) / 45)))) * (V - E_K); % [?A] no capacitance normalization

            % TODO: Incorporate C_m change with stretch??

            current_val = I_Ko;

        end

    end

end