classdef Ins_Buono %< ChannelBase

    % Represents arbitrary TRP channel contributing to NS SAC current in
    % the heart

    % Conducts Na and K
    
    % Does not account for K_o, K_i, Na_o, or Na_i as does original model

    properties (Constant)
        num_gates = 0;
        
        % Values from Buonocunto source code (protocols.m, model_SAC.m)
        g_ns = 10.3203e-6;
        r = 1.4110;
        beta_NS = 0.2203e1; % Differs from Table 1, but matches source code
        ISAC_factor = [4.47, 0.001248]; % Change to [0, 0] for diseased

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
            lambda = X(3); % stretch [%]
            Na_i = X(4);
            K_i = X(5);

            R=8314.0;
            T=310.0;
            F=96485.0;

            nao=140;
            cao=1.8;
            % ko=5.4; %only for I-V curves
            ko=5; %with calcium and contraction
            % ki=142.412524737626; % Initial cond, using as constant
            % nai=12.0996647655188; % Initial cond, using as constant

            ki = K_i;
            nai = Na_i;

            % nai and ki should come from model as they are tracked

            E_Na=(R*T/F)*log(nao/nai);
            E_K=(R*T/F)*log(ko/ki);

            eps = lambda - 1; % strain, maybe rename
            
            % Stretch-dep
            C_m = 1.534e-4 * (((eps^4)/(0.035^4 + eps^4)) + 1); % Output
            % Constant - huge difference, ~0.5x mag as stretch-dep
            % C_m = 1.534e-4;

            sigm_factor = (eps^Ins_Buono.ISAC_factor(1))/(eps^Ins_Buono.ISAC_factor(1)+Ins_Buono.ISAC_factor(2)); % CORRECT
            gamma_SLNs = 0.5 * ((Ins_Buono.beta_NS*eps^3)/(0.15^3 + eps^3));

            I_NsNa = Ins_Buono.r * Ins_Buono.g_ns * sigm_factor * gamma_SLNs * (V-E_Na);
            I_NsK = Ins_Buono.g_ns * sigm_factor * gamma_SLNs * (V-E_K);

            current_val = [I_NsNa, I_NsK];

            % Normalize to capacitance
            current_val = current_val./C_m;

        end

    end

end