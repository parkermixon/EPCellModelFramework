classdef OgiermannParams
    properties(Constant)
        % Markov Transition Parameters
        r1 = 0.02634342; % [1/ms]
        r3 = 0.0008000712; % [1/ms]
        r4 = 0.008945307; % [1/ms]
        r5 = 1.3529515e-5; % [1/ms]
        r6 = 6.206403e-5; % [1/ms]
        r7 = 6.603723e-6; % [1/ms]
        r8 = 0.00079357607; % [1/ms]

        r2 = SAC.OgiermannParams.r1*SAC.OgiermannParams.r3*SAC.OgiermannParams.r5/(SAC.OgiermannParams.r4*SAC.OgiermannParams.r6) % [1/ms]

        ce1 = -2.3753126; % [dimensionless]
        ce3 = -5.6010337; % [dimensionless]
        ce4 = 0.70629007; % [dimensionless]
        cm4 = -12.101605; % [dimensionless]
        ce5 = 6.46389; % [dimensionless]
        ce6 = -1.2173947; % [dimensionless]
        ce7 = 8.739162; % [dimensionless]
        ce8 = -3.364573; % [dimensionless]

        ce2 = SAC.OgiermannParams.ce1 + SAC.OgiermannParams.ce3 + SAC.OgiermannParams.ce5 - SAC.OgiermannParams.ce4 - SAC.OgiermannParams.ce6; % [dimensionless]
        cm2 = -SAC.OgiermannParams.cm4; % [dimensionless]

        % Environment Parameters
        R = 8.314472; % [J/molK]
        T = 308; % [Kelvin]
        F = 96.4853415; % [C/mmol]
        
        K_o = 5.4; % [mM]
        Ca_o = 1.8; % [mM]
        Na_o = 136; % [mM]
        K_i = 140; % [mM]
        prNaK = 0.01833; % [dimensionless] % Free parameter

        % Conductances
        % mahajan_shiferaw_ogiermann_2025.cellml
        scaling = 0.18; % Free parameter
        gPz1K =  0.93; % [uS/nF]
        gPz1Ca = 0.28; % [uS/nF]
        gPz1Na = 1.0; % [uS/nF]
    end
end