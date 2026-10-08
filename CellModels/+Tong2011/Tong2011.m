classdef Tong2011 < CellEPBase

    properties (Constant)

        % V [mV], Ca_i [mM]
        init_state_values = [-53.90915441282156, 0.0001161881607214449];
        
        num_stims = 1; % [Istim]

        channel_classes = struct( ...
            "I_CaL", Tong2011.ICaL, ...
            "I_CaT", Tong2011.ICaT, ...
            "I_Na", Tong2011.INa, ...
            "I_h", Tong2011.Ih, ...
            "I_K1", Tong2011.IK1, ...
            "I_K2", Tong2011.IK2, ...
            "I_Ka", Tong2011.IKa, ...
            "I_KCa", Tong2011.IKCa, ...
            "I_b", Tong2011.Ib, ...
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
        INSCa_index = 11;

        % Important flux indices
        JNaCa_index = 2;
    end

    methods

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