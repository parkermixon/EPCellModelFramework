% Implemented from CellML

classdef IPz1_Ogi < ChannelBase

    properties (Constant)
        num_gates = 4;
        gate_names = {"Pz1_Open", "Pz1_Closed", "Pz1_I1", "Pz1_I2"};

        % r1 = 0.02634342; % [1/ms]
        % r3 = 0.0008000712; % [1/ms]
        % r4 = 0.008945307; % [1/ms]
        % r5 = 1.3529515e-5; % [1/ms]
        % r6 = 6.206403e-5; % [1/ms]
        % r7 = 6.603723e-6; % [1/ms]
        % r8 = 0.00079357607; % [1/ms]
        % 
        % r2 = IPz1_Ogi.r1*IPz1_Ogi.r3*IPz1_Ogi.r5/(IPz1_Ogi.r4*IPz1_Ogi.r6) % [1/ms]
        % 
        % ce1 = -2.3753126; % [dimensionless]
        % ce3 = -5.6010337; % [dimensionless]
        % ce4 = 0.70629007; % [dimensionless]
        % cm4 = -12.101605; % [dimensionless]
        % ce5 = 6.46389; % [dimensionless]
        % ce6 = -1.2173947; % [dimensionless]
        % ce7 = 8.739162; % [dimensionless]
        % ce8 = -3.364573; % [dimensionless]
        % 
        % ce2 = IPz1_Ogi.ce1 + IPz1_Ogi.ce3 + IPz1_Ogi.ce5 - IPz1_Ogi.ce4 - IPz1_Ogi.ce6; % [dimensionless]
        % cm2 = -IPz1_Ogi.cm4; % [dimensionless]
    end

    methods (Static)
        
        function Xg = Init()

            % Xg(1) = 0; % P_Open
            % Xg(2) = 1; % P_Closed
            % Xg(3) = 0; % P_I1
            % Xg(4) = 0; % P_I2

            % ogiermann_2025_piezo1.cellml
            Xg(1) = 0; % P_Open
            Xg(2) = 0.3; % P_Closed
            Xg(3) = 0.1; % P_I1
            Xg(4) = 0.6; % P_I2
            
        end

        function dXg = GetG(G, X, Xg, stim)
    
            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_submem = X(2); % [mM]
            Na_i = X(3); % [mM]
            % TODO: Add K_i instead of getting K_i from Na_i with G.prNaK
            
            % Istim = stim(1);
            p = stim(2); % pressure [mmHg]

            P_Open = Xg(1);
            P_Closed = Xg(2);
            P_I1 = Xg(3);
            P_I2 = Xg(4);

            FonRT = G.F/(G.R*G.T);
            
            gPz1All = G.gPz1K+G.gPz1Ca+G.gPz1Na;
            PPz1Na = G.gPz1Na/gPz1All;
            PPz1Ca = G.gPz1Ca/gPz1All;
            PPz1K   = G.gPz1K/gPz1All;

            eghkPz1 = 1/FonRT*log((PPz1K*G.K_o+PPz1Ca*G.Ca_o+PPz1Na*G.Na_o)/(PPz1K*G.K_i+PPz1Ca*Ca_submem+PPz1Na*Na_i));
            
            gradmu = V - eghkPz1; % [mV]

            % 4 state Markov model
            p_O_I1 = G.r1*exp(G.ce1*gradmu/140);
            p_I1_O = G.r2*exp(G.ce2*gradmu/140+G.cm2*p/70);
            
            p_I1_C = G.r3*exp(G.ce3*gradmu/140);
            p_C_I1 = G.r4*exp(G.ce4*gradmu/140+G.cm4*p/70);

            p_C_O = G.r5*exp(G.ce5*gradmu/140);
            p_O_C = G.r6*exp(G.ce6*gradmu/140);

            p_I2_O = G.r7*exp(G.ce7*gradmu/140*p/70);
            p_O_I2 = G.r8*exp(G.ce8*gradmu/140);

            dP_Open = -(p_O_I1+p_O_I2+p_O_C)*P_Open+p_C_O*P_Closed+p_I1_O*P_I1+p_I2_O*P_I2;
            dP_Closed = -(p_C_I1+p_C_O)*P_Closed+p_O_C*P_Open+p_I1_C*P_I1;
            dP_I1 = -(p_I1_O+p_I1_C)*P_I1+p_C_I1*P_Closed+p_O_I1*P_Open;
            dP_I2 = -p_I2_O*P_I2+p_O_I2*P_Open;

            dXg = zeros(size(Xg));
            dXg(1) = dP_Open;
            dXg(2) = dP_Closed;
            dXg(3) = dP_I1;
            dXg(4) = dP_I2;
        end

        function current_val = GetF(G, X, Xg, stim)

            % Make local copies of state and gating variables
            V = X(1); % [mV]
            Ca_submem = X(2); % [mM]
            Na_i = X(3); % [mM]

            % TODO: Add K_i instead of getting K_i from Na_i with G.prNaK

            p = stim(2); % pressure [mmHg]

            P_Open = Xg(1);
            P_Closed = Xg(2);
            P_I1 = Xg(3);
            P_I2 = Xg(4);         

            FonRT = G.F/(G.R*G.T);

            eks = 1/FonRT*log((G.K_o+G.prNaK*G.Na_o)/(G.K_i+G.prNaK*Na_i));
            ena = 1/FonRT*log(G.Na_o/Na_i);
            ecas = 1/(FonRT*2)*log(G.Ca_o/Ca_submem);

            xiCaPz1 = G.scaling*P_Open*G.gPz1Ca*(V-ecas); % [nA/nF]
            xiNaPz1 = G.scaling*P_Open*G.gPz1Na*(V-ena); % [nA/nF]
            xiKPz1 = G.scaling*P_Open*G.gPz1K*(V-eks); % [nA/nF]
            xiPz1 = xiCaPz1+xiNaPz1+xiKPz1; % [nA/nF]

            current_val = [xiPz1, xiCaPz1, xiNaPz1, xiKPz1];
            
        end

    end

end