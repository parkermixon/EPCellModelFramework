classdef Tong2014 < CellEPBase

    properties (Constant)

        % V [mV], Ca_i [mM]
        init_state_values = [-50.8077440348614, 0.000123555535400408];

        num_stims = 1; % [Istim]

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
            "J_NaCa", Tong2011.JNaCa ...
        );

        % Extra values will be appended to end after channels and fluxes
        % Should be derived from computed currents and fluxes (with no gating variables of their own)
        extra_value_names = {"I_NaCa"}

        % Important channel indices
        ICaL_index = 1;
        ICaT_index = 2;
        INSCa_index = 14;

        % Important flux indices
        JNaCa_index = 2;
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

            values = I_NaCa_value;
        end

        function dX = ComputeStateDerivatives(obj, G, X, stim, channel_values, flux_values, extra_values)

            dX = zeros(size(X));

            % Handle J_Camem
            I_CaL_value = channel_values(obj.ICaL_index);
            I_CaT_value = channel_values(obj.ICaT_index);
            I_NSCa_value = channel_values(obj.INSCa_index);

            % Compute derivatives
            % dV/dt
            dX(1) = -(sum(channel_values) + sum(extra_values) + stim);

            % dCa_i/dt
            J_Camem = ((G.AV_c*G.C_m*G.bet)/(G.z_Ca*G.F))*(I_CaL_value+I_CaT_value+I_NSCa_value);

            dX(2) = -(J_Camem + sum(flux_values));

        end

    end

end