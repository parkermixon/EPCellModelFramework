classdef Tong2011_SAC < CellEPBase

    properties (Constant)

        state_names = {"V", "Ca_i"};
        % V [mV], Ca_i [mM]
        init_state_values = [-53.90915441282156, 0.0001161881607214449];

        num_stims = 2; % [Istim, lambda]

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
            "I_NaK", Tong2011.INaK, ...
            "I_Pz1", Tong2011_SAC.IPz1_Wrapper, ... % Assumed to only be Ca2+, although is actually NS
            "I_TK1", Tong2011_SAC.ITK1_Wrapper... % K+ SAC
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
        IPz1_index = 16;

        % Important flux indices
        JNaCa_index = 2;

    end

    methods % Copied from Tong2011

        function values = ComputeExtraValues(obj, G, X, channel_values, flux_values)
            JNaCa_value = flux_values(obj.JNaCa_index);
            I_NaCa_value = -0.5 * ((G.z_Ca*G.F*JNaCa_value)/(G.C_m*G.bet)) * (1/G.AV_c);

            values = I_NaCa_value;
        end

        % Modified to include Ca2+ flux through IPz1
        function dX = ComputeStateDerivatives(obj, G, X, stim, channel_values, flux_values, extra_values)

            dX = zeros(size(X));

            Istim = stim(1);

            % Handle J_Camem
            I_CaL_value = channel_values(obj.ICaL_index);
            I_CaT_value = channel_values(obj.ICaT_index);
            I_NSCa_value = channel_values(obj.INSCa_index);
            I_Pz1_value = channel_values(obj.IPz1_index);

            % Compute derivatives
            % dV/dt
            dX(1) = -(sum(channel_values) + sum(extra_values) + Istim);

            % dCa_i/dt
            J_Camem = ((G.AV_c*G.C_m*G.bet)/(G.z_Ca*G.F))*(I_CaL_value+I_CaT_value+I_NSCa_value+I_Pz1_value);

            dX(2) = -(J_Camem + sum(flux_values));

        end

        function stim_mod = ConvertStimulus(stim)

            %% TODO: What type of pressure is this? Balloon pressure or pressure clamp
            % And what type of pressure is used in other Piezo1 models that they mention in their papers??
            % Do we care?? They are already converting stretch into the stimulus used in their model

            % Comment from OpenCOR Ogiermann 2026 model "mixed_pacing_experiment_QuinnKohl2016.cell"
            % Conversion of stretch to pressure with some guesstimated scaling factor. 
            % This likely does not reflect reality and needs further investigation in future iterations of this model.
            % Hence, we could also just apply the pressure directly here.
            % If we assume that a stretch of 30% saturates the ion channel 
            % response (70 mmHg in the experiment) and a linear relation between stretch and pressure, we get the factor of 350 present here.

            Istim = stim(1);
            lambda = stim(2);
            
            % Convert lambda to pressure following Ogiermann et al. 2026 (OpenCOR)
            % Modified 250->350 to match comment???
            p = 350.0*max(0, lambda-1); % [N or mmHg]??
            
            stim_mod = [Istim, lambda, p];

        end

        % Wrap public methods that pass stim to convert lambda stimulus input to pressure
        function [O, Ga, D] = Run(obj, G, stim, t_range, isVStim)

            stim_mod = obj.ConvertStimulus(stim);

            [O, Ga, D] = Run@CellEPBase(obj, G, stim_mod, t_range, isVStim);

        end

        function [O, Ga, D] = RunODE15s(obj, G, stim_func, tspan, isVStim, rel_tol, abs_tol, max_step)
            
            stim_func_mod = @(t) obj.ConvertStimulus(stim_func(t));

            [O, Ga, D] = RunODE15s@CellEPBase(obj, G, stim_func_mod, tspan, isVStim, rel_tol, abs_tol, max_step);

        end

    end
end