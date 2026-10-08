classdef Yang2026 < CellEPBase

    % Implementing OXT -> IP3 -> Tong 2014 portions

    properties (Constant)

        %TODO: C_oxt should be stim

        % V [mV], Ca_i [mM], C_oxt [nM], C_oxtr [nM], C_oxtrc [nM], C_oxtrc_int [nM], IP3 [nM]
        state_names = {"V", "Ca_i", "IP3", "C_oxtr", "C_oxtrc", "C_oxtrc_int"};
        % init_state_values = [-49, 0.0002, 43.59, 216, 0.0111, 0.0272]; % src, ish?

        % C_oxtr value testing, Steady State values after 100s no stim
        % (close enough)
        % C_oxtr = 1 nM
        % init_state_values = [-51.3670, 1.0917e-04, 0.0894, 1, 1.2005e-07, 0.0011];
        % C_oxtr = 10 nM
        % init_state_values = [-51.3670, 1.0917e-04, 2.5272, 10, 3.6503e-05, 0.0060];
        % C_oxtr = 100 nM
        % init_state_values = [-51.3670, 1.0917e-04, 0.0894, 100, 1.2005e-07, 0.0011];
        % C_oxtr = 200 nM
        init_state_values = [-51.3670, 1.0917e-04, 0.0030, 200, 4.0589e-10, 2.1290e-04];

        % init_state_values = [-53.222896781073295, 1.940418038357952e-05] % evaluated @ IP3 = 0
        
        num_stims = 2; % [Istim, C_oxt]

        channel_classes = struct( ...
            "I_CaL", Tong2011.ICaL, ...
            "I_CaT", Tong2011.ICaT, ...
            "I_Na", Tong2011.INa, ...
            "I_h", Tong2011.Ih, ...
            "I_K1", Tong2011.IK1, ...
            "I_KCNQ1", Tong2014.IKCNQ1, ...
            "I_KCNQ4", Tong2014.IKCNQ4, ...
            "I_KCNQ5", Tong2014.IKCNQ5, ...
            "I_hERG", Tong2014.IhERG, ...
            "I_K2", Tong2011.IK2, ...
            "I_Ka", Tong2011.IKa, ...
            "I_KCa", Tong2011.IKCa, ...
            "I_ClCa", Tong2011.IClCa, ...
            "I_NSCa", Tong2011.INSCa, ...
            "I_NSK", Tong2011.INSK, ...
            "I_NSL", Tong2011.INSL, ...
            "I_NSNa", Tong2011.INSNa, ...
            "I_NaK", Tong2011.INaK ...
        );

        flux_classes = struct( ...
            "J_PMCA", Tong2011.JPMCA, ...
            "J_NaCa", Tong2011.JNaCa, ...
            "J_rel", Yang2026.Jrel, ... % Yang2026.Jrel_Gated
            "J_up", Yang2026.Jup, ...
            "J_SOC", Yang2026.JSOC ...
        );

        % Extra values will be appended to end after channels and fluxes
        % Should be derived from computed currents and fluxes (with no gating variables of their own)
        extra_value_names = {"I_NaCa", "I_eq", "J_t"}

        % Important channel indices
        ICaL_index = 1;
        ICaT_index = 2;
        INSCa_index = 14;

        % Important flux indices
        JNaCa_index = 2;
        Jrel_index = 3;
        Jup_index = 4;
        JSOC_index = 5;
    end

    methods

        function [X, Xg] = Init(obj, init_values)
            arguments
                obj CellEPBase
                init_values struct = struct()
            end

            default_init_values = struct();

            % Update default Tong 2011 values to Tong 2014 values (from OpenCOR)
            % I_CaL, [d, f_1, f_2]
            default_init_values.I_CaL = [0.0160574909192411, 0.861723329545759, 0.861723329545101];
            % I_CaT, [b, g]
            default_init_values.I_CaT = [0.585739993599288, 0.0281751834106473];
            % I_ClCa, [c]
            default_init_values.I_ClCa = 0.00066950904540682;
            % I_h, [y]
            default_init_values.I_h = 0.00182158784678185;
            % I_K1, [q, r_1, r_2]
            default_init_values.I_K1 = [0.22560249352574, 0.126167455396881, 0.126167442021144];
            % I_K2, [p, k_1, k_2]
            default_init_values.I_K2 = [0.135874773287517, 0.994482738453784, 0.994482738453784];
            % I_Ka, [s, x]
            default_init_values.I_Ka = [0.0456227251358283, 0.0616278982343972];
            % I_KCa, [x_a, x_ab1]
            default_init_values.I_KCa = [0.00035751229735921, 0.00267392787579562];
            % I_Na, [m, h]
            default_init_values.I_Na = [0.166998688814229, 0.315607550727852];

            % Copy over supplied init values
            if ~isempty(init_values)
                for f = fieldnames(init_values)'
                    default_init_values.(f) = init_values.(f);
                end
            end

            [X, Xg] = Init@CellEPBase(obj, default_init_values);

        end

        function values = ComputeExtraValues(obj, G, X, channel_values, flux_values)
            JNaCa_value = flux_values(obj.JNaCa_index);
            I_NaCa_value = -0.5 * ((G.z_Ca*G.F*JNaCa_value)/(G.C_m*G.bet)) * (1/G.AV_c);

            % SR dynamics
            Jrel_value = flux_values(obj.Jrel_index);
            Jup_value = flux_values(obj.Jup_index);
            JSOC_value = flux_values(obj.JSOC_index);

            J_t_value = Jrel_value + JSOC_value - Jup_value;

            I_eq = -G.z_Ca * G.F * J_t_value / (G.AV_c * G.C_m * G.bet);

            values = [I_NaCa_value, I_eq, J_t_value];
        end

        function dX = ComputeStateDerivatives(obj, G, X, stim, channel_values, flux_values, extra_values)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]
            IP3 = X(3); % [nM]
            C_oxtr = X(4); % [nM]
            C_oxtrc = X(5); % [nM]
            C_oxtrc_int = X(6); % [nM]
            
            Istim = stim(1);
            C_oxt = stim(2);

            % Initialize dX to 0
            dX = zeros(size(X));

            % Handle J_Camem
            I_CaL_value = channel_values(obj.ICaL_index);
            I_CaT_value = channel_values(obj.ICaT_index);
            I_NSCa_value = channel_values(obj.INSCa_index);

            % Compute derivatives
            % dV/dt
            dX(1) = -(sum(channel_values) + sum(extra_values(1:2)) + Istim);

            % dCa_i/dt
            J_Camem = ((G.AV_c*G.C_m*G.bet)/(G.z_Ca*G.F))*(I_CaL_value+I_CaT_value+I_NSCa_value);
            
            % Paper does not state if SR fluxes are counted in dCa_i, we are
            % Jrel_value = flux_values(obj.Jrel_index);
            % Jup_value = flux_values(obj.Jup_index);
            % JSOC_value = flux_values(obj.JSOC_index);
            % 
            % J_t_value = Jrel_value + JSOC_value - Jup_value;

            dX(2) = -(J_Camem + sum(flux_values(1:2))) + extra_values(3);

            % dIP3/dt
            Ca_i_nM = Ca_i*1e6; % [mM] -> [nM]

            % These are intermediate values include in the full IP3 eq (coeff1)
            PLCbeta_act = G.PLCbeta_tot * (C_oxtrc / (C_oxtrc + G.K_d_PLCbeta));
            v_PLC = G.k_cat * PLCbeta_act;

            % coeff1 = (G.k_cat * C_oxtrc) / (C_oxtrc + G.K_d_PLCbeta);
            coeff2 = G.v_5 / (1 + G.k_5_IP3/IP3);
            coeff3 = G.v_6 / ((1 + G.k_6_IP3/IP3) * (1 + G.k_7_IP3/Ca_i_nM));

            % Text states it as direct value, Yang & Jafri instead have as ODE
            dIP3_dt = v_PLC - coeff2 - coeff3;
            % X(3) = IP3_n;
            dX(3) = dIP3_dt;

            % OXTR dynamics
            % Evaluate rate constants
            k_int = (G.k_int_max * C_oxtrc) / (G.k_phos + C_oxtrc);
            k_rec = k_int;

            % OXTR Pathways
            OXTR_OXTRC_on = G.k_on * C_oxt * C_oxtr;
            OXTR_OXTRC_off = G.k_off * C_oxtrc;
            OXTR_OXTRC_INT_rec = k_rec * C_oxtrc_int;
            OXTR_OXTRC_INT_int = k_int * C_oxtrc;

            dCoxtr_dt = -OXTR_OXTRC_on + OXTR_OXTRC_off + OXTR_OXTRC_INT_rec;
            dX(4) = dCoxtr_dt;

            dCoxtrc_dt = OXTR_OXTRC_on - OXTR_OXTRC_off - OXTR_OXTRC_INT_int - (G.k_des * C_oxtrc);
            dX(5) = dCoxtrc_dt;

            dCoxtrc_int_dt = OXTR_OXTRC_INT_int - OXTR_OXTRC_INT_rec;
            dX(6) = dCoxtrc_int_dt;  
        end

    end

end