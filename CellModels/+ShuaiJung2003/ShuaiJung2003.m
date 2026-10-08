classdef ShuaiJung2003 < CellEPBase

    properties (Constant)

        % Ca_i [uM]
        init_state_values = [0.017575686927209];
        
        num_stims = 1; % [IP3]

        channel_classes = struct();

        flux_classes = struct( ...
            "J_Pump", ShuaiJung2003.JPump, ...
            "J_Leak", ShuaiJung2003.JLeak, ...
            "J_Channel", ShuaiJung2003.JChannel ...
        );

        J_Pump_idx = 1;
        J_Leak_idx = 2;
        J_Channel_idx = 3;

        extra_value_names = {};
    end

    methods

        function values = ComputeExtraValues(obj, G, X, channel_values, flux_values)
            values = [];
        end

        function dX = ComputeStateDerivatives(obj, G, X, stim, channel_values, flux_values, extra_values)

            dX = zeros(size(X));

            J_Pump = flux_values(obj.J_Pump_idx);
            J_Leak = flux_values(obj.J_Leak_idx);
            J_Channel = flux_values(obj.J_Channel_idx);

            % Compute derivatives
            % dCa_i/dt
            dX(1) = -J_Pump + J_Leak + J_Channel;

        end

    end

end