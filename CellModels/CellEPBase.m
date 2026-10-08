classdef CellEPBase

    properties (Abstract, Constant)
        % Cell Model specific variables
        state_names (1, :) cell;
        init_state_values (1, :) double;
        channel_classes struct;
        flux_classes struct;
        % Extra values will be appended to end after channels and fluxes
        % Should be derived from computed currents and fluxes (with no gating variables of their own)
        extra_value_names (1, :) cell;
        num_stims (1, 1) double; % Istim always 1st index
    end

    properties
        % Helper variables/constants, computed from variables above
        num_state_vars (1, 1) double;
        channel_names (:, 1) cell;
        channel_gate_counts (:, 1) double;
        
        flux_names (:, 1) cell;
        flux_gate_counts (:, 1) double;
        
        total_gate_counts (:, 1) double;

        num_extra_values (1, 1) double;
    end

    methods (Abstract)
        
        values = ComputeExtraValues(obj, G, X, channel_values, flux_values);

        dX = ComputeStateDerivatives(obj, G, X, stim, channel_values, flux_values, extra_values);

    end

    methods

        function obj = CellEPBase()

            % Ensure num_state_vars equals number initial state values
            assert(length(obj.state_names) == length(obj.init_state_values), "State array init lengths don't match");

            % Compute and set helper variables
            obj.num_state_vars = length(obj.init_state_values);

            obj.channel_names = fieldnames(obj.channel_classes);
            obj.channel_gate_counts = cellfun(@(c) c.num_gates, struct2cell(obj.channel_classes));

            obj.flux_names = fieldnames(obj.flux_classes);
            obj.flux_gate_counts = cellfun(@(c) c.num_gates, struct2cell(obj.flux_classes));

            obj.total_gate_counts = vertcat(obj.channel_gate_counts, obj.flux_gate_counts);

            obj.num_extra_values = length(obj.extra_value_names);

        end

        function [X, Xg] = Init(obj, init_values)
            arguments
                obj CellEPBase
                init_values struct = struct()
            end

            X = zeros(1, obj.num_state_vars);
            X(:) = obj.init_state_values;

            % Setup Xg as a cell array to more easily associate gating
            % variables with each current
            Xg = {};

            % Channels
            for channel_name = string(fieldnames(obj.channel_classes))'
                fprintf("Initializing %s\n", channel_name);
                
                % Default values from class
                Xg_channel = obj.channel_classes.(channel_name).Init();

                if isfield(init_values, channel_name)
                    
                    disp("Found given values")

                    if length(init_values.(channel_name)) == obj.channel_classes.(channel_name).num_gates
                        Xg_channel = init_values.(channel_name);
                    else
                        disp("Error, mismatch in number of gates. Using default values.")
                    end
                end
                
                disp(Xg_channel)
                Xg = vertcat(Xg, {Xg_channel});
            end

            % Fluxes
            for flux_name = string(fieldnames(obj.flux_classes))'
                fprintf("Initializing %s\n", flux_name);

                % Default values from class
                Xg_flux = obj.flux_classes.(flux_name).Init();

                if isfield(init_values, flux_name)

                    disp("Found given values")

                    if length(init_values.(flux_name)) == obj.flux_classes.(flux_name).num_gates
                        Xg_flux = init_values.(flux_name);
                    else
                        disp("Error, mismatch in number of gates. Using default values.")
                    end
                end

                disp(Xg_flux)
                Xg = vertcat(Xg, {Xg_flux});
            end
        end

        function [O, Ga, D] = Run(obj, G, stim, t_range, isVStim)
            % NOTE: does not return initial current values to output
            arguments
                obj CellEPBase
                G (1, 1) struct
                stim (:, :) double % col count should match t_range length
                t_range (1, :) double
                isVStim (1, 1) = false
            end

            dt = t_range(2)-t_range(1);
            N_total = length(t_range);

            % Check that stim array is correct size
            assert(size(stim, 1) == obj.num_stims, "Number of stimuli must match num_stims");
            assert(size(stim, 2) == N_total, "Number stim cols ~= t_range length");

            % Create and initialize X and Xg arrays

            [X, Xg] = obj.Init();

            values = zeros(1, length(obj.channel_names)...
                + length(obj.flux_names) ...
                + obj.num_extra_values ...
                );

            % Create output arrays
            outputs = zeros(N_total, length(X)+length(values));
            outputs(1, :) = horzcat(X, values);

            gate_outputs = zeros(N_total, sum(obj.total_gate_counts));
            gate_outputs(1, :) = [Xg{:}];

            deltas = zeros(N_total, length(X) + sum(obj.total_gate_counts));
            deltas(1, :) = 0;

            for i=1:N_total
                if (isVStim)
                    X(1) = stim(i);

                    [dX, dXg, values] = obj.Update(G, X, Xg, 0);

                    dX(1) = 0;
                else
                    [dX, dXg, values] = obj.Update(G, X, Xg, stim(:, i));
                end

                [X, Xg] = obj.ForwardEulerODE(X, Xg, dX, dXg, dt);

                outputs(i, :) = horzcat(X, values);
                gate_outputs(i, :) = [Xg{:}];
                % deltas(i, :) = horzcat(dX, dXg); % Not yet implemented consistently
            end

            [O, Ga, D] = obj.CreateOutputStructs(outputs, gate_outputs, deltas, t_range);
        end
    
        function [O, Ga, D] = RunODE15s(obj, G, stim_func, tspan, isVStim, rel_tol, abs_tol, max_step)
            arguments
                obj CellEPBase
                G (1, 1) struct
                stim_func (1, 1) 
                tspan (1, 2) double
                isVStim (1, 1) = false
                rel_tol (1, 1) double = 1e-7
                abs_tol (1, 1) double = 1e-9
                max_step (1, 1) double = 10
            end

            % Create and initialize X and Xg arrays

            [X, Xg] = obj.Init();

            % Linearize Xg
            Xg_lin = [Xg{:}];

            % Run using Matlab ode solver
            y0 = horzcat(X, Xg_lin).';

            opts = odeset('RelTol', rel_tol,'AbsTol', abs_tol,'MaxStep', max_step);

            X_idxs = 1:length(X);
            Xg_idxs = (length(X)+1):(length(X)+length(Xg_lin));
            
            [t, y] = ode15s(@(t, y) obj.UpdateODEWrapper(t, y, G, stim_func, isVStim, X_idxs, Xg_idxs), tspan, y0, opts);

            % Prepare outputs
            X_vals = y(:, 1:obj.num_state_vars);
            Xg_vals_lin = y(:, obj.num_state_vars+1:end);

            % Convert Xg to cell array (should be moved to a sep util function)
            Xg_vals_cell = mat2cell(Xg_vals_lin, [length(t)], obj.total_gate_counts);

            Xg_offset = 0;

            % Retrieve current values - TODO: clean up, make more efficient
            % possibly (don't need so many loops over t, can be combined)
            % TODO: Move this to CreateOutputStructs, to be done if no channel values provided
            % THIS IS VERY SIM TO UPDATE AND CREATEOUTPUT...
            channel_values = zeros(length(t), length(obj.channel_names));

            for i = 1:length(obj.channel_names)

                channel_class = obj.channel_classes.(obj.channel_names{i});

                Xg_channel = Xg_vals_cell{Xg_offset + i};

                for n = 1:length(t)
                    channel_values(n, i) = channel_class.GetF(G, X_vals(n, :), Xg_channel(n, :), stim_func(t(n)));
                end
            end

            % Update offset to account for channels
            Xg_offset = Xg_offset + length(obj.channel_names);

            % Retrieve flux values
            flux_values = zeros(length(t), length(obj.flux_names));

            for i = 1:length(obj.flux_names)

                flux_class = obj.flux_classes.(obj.flux_names{i});

                Xg_flux = Xg_vals_cell{Xg_offset + i};

                for n = 1:length(t)
                    flux_values(n, i) = flux_class.GetF(G, X_vals(n, :), Xg_flux(n, :), stim_func(t(n)));
                end
            end

            % Retrieve extra values, should have no gates
            extra_values = zeros(length(t), obj.num_extra_values);

            for n = 1:length(t)

                extra_values(n, :) = obj.ComputeExtraValues(G, X_vals(n, :), channel_values(n, :), flux_values(n, :));

            end

            % Create output structs
            outputs = horzcat(X_vals, channel_values, flux_values, extra_values); % extra values appended to end
            gate_outputs = Xg_vals_lin;
            deltas = []; % deltas is not implemented, also not really worth it to implement

            [O, Ga, D] = obj.CreateOutputStructs(outputs, gate_outputs, deltas, t);
        end
    
        function [dX, dXg, values] = Update(obj, G, X, Xg, stim)

            % Initialize output arrays
            dX = zeros(size(X));
            dXg = cell(size(Xg));

            Xg_offset = 0;

            % Update channels
            channel_values = zeros(1, length(obj.channel_names));

            for i = 1:length(obj.channel_names)
                channel_class = obj.channel_classes.(obj.channel_names{i});
                Xg_channel = Xg{i};

                % Get current current value
                value = channel_class.GetF(G, X, Xg_channel, stim);

                % Get gating variable derivatives based on current current values
                dXg_channel = channel_class.GetG(G, X, Xg_channel, stim);

                dXg{i} = dXg_channel;
                channel_values(i) = value;
            end

            % Update offset to account for channels
            Xg_offset = Xg_offset + length(obj.channel_names);

            % Update fluxes
            flux_values = zeros(1, length(obj.flux_names));

            for i = 1:length(obj.flux_names)
                flux_class = obj.flux_classes.(obj.flux_names{i});
                Xg_flux = Xg{Xg_offset + i};

                % Get current current value
                value = flux_class.GetF(G, X, Xg_flux, stim);

                % Get gating variable derivatives based on current current values
                dXg_flux = flux_class.GetG(G, X, Xg_flux, stim);

                dXg{Xg_offset + i} = dXg_flux;
                flux_values(i) = value;
            end

            % Compute extra values
            extra_values = obj.ComputeExtraValues(G, X, channel_values, flux_values);

            % Compute state variable time derivatives
            dX = obj.ComputeStateDerivatives(G, X, stim, channel_values, flux_values, extra_values);

            % List values to output
            % Extra values are appended to end
            values = horzcat(channel_values, flux_values, extra_values);
        end

        function dy = UpdateODEWrapper(obj, t, y, G, stim_func, isVStim, X_idxs, Xg_idxs)

            X = y(X_idxs);
            % Convert Xg to cell array used by Update
            Xg_lin = y(Xg_idxs);
            Xg = mat2cell(Xg_lin, obj.total_gate_counts, [1]);

            if (isVStim)
                X(1) = stim_func(t);
                % Set Istim to 0 but retain other stim values
                stim_vals = stim_func(t);
                stim_vals(1) = 0;
                [dX, dXg, ~] = obj.Update(G, X, Xg, stim_vals);
                dX(1) = 0;
            else
                [dX, dXg, ~] = obj.Update(G, X, Xg, stim_func(t));
            end

            dy = zeros(size(y));
            dy(X_idxs) = dX;
            dy(Xg_idxs) = [dXg{:}]; % Linearize dXg

        end

        function [O, Ga, D] = CreateOutputStructs(obj, outputs, gate_outputs, deltas, t_range)

            O = struct();
            Ga = struct();
            D = struct(); % Not currently supported

            O.Time = t_range;

            % O.V = outputs(:, 1);
            % O.Ca_i = outputs(:, 2);

            % Output state values
            for i = 1:obj.num_state_vars
                O.(obj.state_names{i}) = outputs(:, i);
            end
            

            gate_outputs_cell = mat2cell(gate_outputs, [length(t_range)], obj.total_gate_counts);

            % Much copied from Update(), should maybe be moved to a
            % function

            output_offset = obj.num_state_vars;

            % Set channel and respective gate outputs
            for i = 1:length(obj.channel_names)

                % Set channel values
                O.(obj.channel_names{i}) = outputs(:, output_offset+i); 

                channel_class = obj.channel_classes.(obj.channel_names{i});

                Xg_channel = gate_outputs_cell{i};

                for g = 1:size(Xg_channel, 2)
                    Ga.(channel_class.gate_names{g}) = Xg_channel(:, g);
                end
            end

            % Update offset to account for channels
            output_offset = output_offset + length(obj.channel_names);
            gate_offset = length(obj.channel_names);

            % Set fluxes and respective gate outputs
            for i = 1:length(obj.flux_names)

                % Set channel values
                O.(obj.flux_names{i}) = outputs(:, output_offset+i); 

                flux_class = obj.flux_classes.(obj.flux_names{i});

                Xg_flux = gate_outputs_cell{gate_offset+i};

                for g = 1:size(Xg_flux, 2)
                    Ga.(flux_class.gate_names{g}) = Xg_flux(:, g);
                end
            end

            % Update offset to account for fluxes
            output_offset = output_offset + length(obj.flux_names);

            % Set extra outputs; should not have any gating
            for i = 1:obj.num_extra_values
                O.(obj.extra_value_names{i}) = outputs(:, output_offset+i);
            end
        end

    end

    methods (Static)

        function [X, Xg] = ForwardEulerODE(nX, nXg, dX, dXg, dt)
            X = nX + dt*dX;

            % Xg = nXg + dt*dXg; but in cell array form
            Xg = cellfun(@(nc, dc) nc + dt*dc, nXg, dXg, 'UniformOutput', false);
        end

    end

end