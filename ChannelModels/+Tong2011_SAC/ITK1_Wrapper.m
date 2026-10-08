classdef ITK1_Wrapper < SAC.IKo_Buono

    % Wrapper for Buonocunto et al. TREK-1 (IKo) model

    methods (Static)

        % Update gating variables
        function dXg = GetG(G, X, Xg, stim)

            % Convert inputs to form used by IKo_Buono
            [X_mod, stim_mod] = Tong2011_SAC.ITK1_Wrapper.WrapInputs(G, X, stim);

            % Pass to IKo_Buono
            dXg = GetG@SAC.IKo_Buono(G, X_mod, Xg, stim_mod);
        end

        % Update current value
        function current_val = GetF(G, X, Xg, stim)

            % Convert inputs to form used by IKo_Buono
            [X_mod, stim_mod] = Tong2011_SAC.ITK1_Wrapper.WrapInputs(G, X, stim);

            % Pass to IKo_Buono
            ITK1 = GetF@SAC.IKo_Buono(G, X_mod, Xg, stim_mod);

            % Return current value
            current_val = ITK1;
        end

        % Wrap inputs to change to form used by IPz1_Ogi
        function [X_mod, stim_mod] = WrapInputs(G, X, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_i = X(2); % [mM]

            % Pressure is passed along with lambda
            Istim = stim(1); % [pA/pF]
            lambda = stim(2);
            % p = stim(3);

            X_mod = [V, Ca_i];
            stim_mod = [Istim, lambda];

        end
    end
end