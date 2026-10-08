classdef IPz1_Wrapper < SAC.IPz1_Ogi

    % Wrapper for Ogiermann et al. Piezo1 model

    methods (Static)

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Convert inputs to form used by IPz1_Ogi
            [X_mod, stim_mod] = Tong2011_SAC.IPz1_Wrapper.WrapInputs(G, X, stim);

            % Pass to IPz1_Ogi
            dXg = GetG@SAC.IPz1_Ogi(G, X_mod, Xg, stim_mod);
        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Convert inputs to form used by IPz1_Ogi
            [X_mod, stim_mod] = Tong2011_SAC.IPz1_Wrapper.WrapInputs(G, X, stim);

            % Pass to IPz1_Ogi
            IPz1 = GetF@SAC.IPz1_Ogi(G, X_mod, Xg, stim_mod);

            % Only returning total current
            current_val = IPz1(1);
        end

        % Wrap inputs to change to form used by IPz1_Ogi
        function [X_mod, stim_mod] = WrapInputs(G, X, stim)
            
            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            % Pressure is passed along with lambda
            Istim = stim(1);
            % lambda = stim(2);
            p = stim(3);

            % Convert to order used by IPz1
            Na_i = G.Na_i; % [mM]
            Ca_submem = Ca_i; % We place Piezo-1 in main cytosol

            X_mod = [V, Ca_submem, Na_i];
            stim_mod = [Istim, p];

        end
    end
end