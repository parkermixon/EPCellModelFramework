classdef BuonocuntoParams
    properties(Constant)
        g_Ko = 12.8208e-6;
        beta_Ko = 0.8742;

        R=8314.0;
        T=310.0;
        F=96485.0;
        
        % Environment Parameters
        ko=5.4; % Used in their model only for I-V curves
        ki=142.412524737626; % Their initial cond, using as constant
    end
end