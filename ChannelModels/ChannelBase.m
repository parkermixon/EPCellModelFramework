classdef ChannelBase

    properties (Abstract, Constant)
        num_gates;
        gate_names;
    end

    methods (Static, Abstract)
        
        Xg = Init();
        
        % Update gating variables
        dXg = GetG(G, X, Xg, stim);

        % Update current value
        current_val = GetF(G, X, Xg, stim);

    end

end