% All values originally from Table S4 in paper; then verified using source
% code and CellML model

classdef OriginalParams
    properties(Constant)
        % T: Temperature
        T = 308; % [K] - CellML
        % F: Faraday constant
        F = 96485; % [C/mol] - CellML
        % R: Gas constant
        R = 8314; % [J/K/kmol] - CellML
        % Na+ valency
        z_Na = 1; % [dimensionless] - CellML
        % K+ valency
        z_K = 1; % [dimensionless] - CellML
        % Ca2+ valency
        z_Ca = 2; % [dimensionless] - CellML
        % Intracellular Na+ concentration
        Na_i = 4; % [mM] - CellML
        % Extracellular Na+ concentration
        Na_o = 130; % [mM] - CellML
        % Intracellular K+ concentration
        K_i = 140; % [mM] - CellML
        % Extracellular K+ concentration
        K_o = 6; % [mM] - CellML
        % Intracellular Cl- concentration
        Cl_i = 46; % [mM] - CellML
        % Extracellular Cl- concentration
        Cl_o = 130; % [mM] - CellML
        % Extracellular Ca+ concentration
        Ca_o = 2.5; % [mM] - CellML
        % Extracellular Mg2+ concentration
        Mg_o = 0.5; % [mM] - CellML
        % Specific membrane capacitance
        C_m = 1; % [uF/cm^2] - CellML
        % Cell surface area to volume ratio
        AV_c = 4; % [1/m] - CellML
        % Proportion of free Ca2+ ions
        % defined as buff in CellML & src code
        bet = 0.015; % [dimensionless] - CellML
        
        % Plasma membrane Ca-ATPase
        % PMCA Ca2+ flux
        J_PMCA = 3.5e-7; % [mM/ms] - CellML
        % Half-saturation concentration for PMCA
        % Defined in uM, converted to mM
        K_mPMCA = 0.0005; % [mM] - CellML
        % Hill coefficient for PMCA
        n_PMCA = 2; % [dimensionless] - CellML
        
        % NaCa Exchanger
        % NaCa exchanger flux
        J_NaCa = 3.5e-6; % [mM/ms] - CellML
        % Half-saturation concentration for NaCa exchanger allosteric factor
        % Wrong in Table S4; corrected with src code and CellML model
        K_mAllo = 0.0003; % [mM] - CellML
        % Hill coefficient for NaCa exchanger allosteric factor
        n_Allo = 4; % [dimensionless] - CellML
        % Saturation factor for NaCa exchanger
        k_sat = 0.27; % [dimensionless] - CellML
        % Partition parameter for NaCa exchanger flux
        % xgamma in CellML
        gam = 0.35; % [dimensionless] - CellML
        % Na_i dissociation constant for NaCa exchanger
        K_mNai = 30.0; % [mM] - CellML
        % Ca_i dissociation constant for NaCa exchanger
        % Defined in uM, converted to mM; checked with CellML & src code
        K_mCai = 7e-3; % [mM] - CellML
        % Na_o dissociation constant for NaCa exchanger
        K_mNao= 87.5; % [mM] - CellML
        % Ca_o dissociation constant for NaCa exchanger
        K_mCao = 1.3; % [mM] - CellML
        
        % L-type Calcium (ICaL)
        % Maximum conductance of ICaL
        g_CaL = 0.6; % [mS/uF] - CellML
        % Reversal potential of ICaL
        E_CaL = 45; % [mV] - CellML
        % Half-saturation concentration for ICaL
        % Defined in uM, converted to mM
        % kmca in CellML
        K_mCaL = 1e-3; % [mM] - CellML
        % Inactivation Steady State Slope Factor
        f_ss_sf = 7; % [mV] - CellML / src code
        
        % T-type Calcium (ICaT)
        % Maximum conductance of ICaT
        g_CaT = 0.058; % [mS/uF] - CellML
        % Reversal potential of ICaT
        E_CaT = 42; % [mV] - CellML
        
        % Sodium Current (INa)
        % maximum conductance of INa
        % can range 0 to 0.12 (double check)
        g_Na = 0; % [mS/uF] - CellML
        
        % Hyperpolarization-activated Current (Ih)
        % maximum conductance of Ih
        g_h = 0.0542; % [mS/uF] - CellML
        % Permeability of I_h to Na:K
        % defined as P_Na=0.35 and P_K=1 separately in CellML; ratio is the same
        P_NaK = 0.35; % [dimensionless]
        
        % Potassium Currents (IK1, IK2, IKa, IKb, IK(Ca))
        % Maximum conductance of total IK
        % Not used in CellML; instead the g_K1, etc. expressions are solved
        g_K = 0.8; % [mS/uF]
        % Maximum conductance of background potassium current
        g_b = 0.004; % [mS/uF] - CellML
        % Maximum conductance of IK1
        g_K1 = 0.24; % [mS/uF] - CellML
        % Maximum conductance of IhERG
        g_hERG = 0.08; % [mS/uF] - CellML
        % Maximum conductance of KCNQ1
        g_KCNQ1 = 0.0032; % [mS/uF] - CellML
        % Maximum conductance of KCNQ4
        g_KCNQ4 = 0.024; % [mS/uF] - CellML
        % Maximum conductance of KCNQ5
        g_KCNQ5 = 0.016; % [mS/uF] - CellML
        % Maximum conductance of IK2
        g_K2 = 0.032; % [mS/uF] - CellML
        % Maximum conductance of IKa
        g_Ka = 0.16; % [mS/uF] - CellML
        % ???
        % g_BK = OriginalParams.g_K; % Never used, but defined in Table S4
        % Maximum conductance of calcium-activated potassium currents (BK, BKaB)
        g_KCa = 0.8; % [mS/uF] - CellML
        % proportion of I_a in IK(Ca)
        % Value from paper source code; not defined in Table S4
        % gbka in CellML
        p_a = 0.2; % [dimensionless] - CellML
        % proportion of I_ab1 in IK(Ca)
        % Value from paper source code; not defined in Table S4
        % gbkab in CellML
        p_b = 0.1;  % [dimensionless] - CellML
        
        % Nonselective Cation Current (I_NSCC)
        % Maximum conductance of I_NSCC leak current
        % Used here but not defined in CellML
        g_L = 0; % [mS/uF] % 0.009685 mS/pF by Testrow et al. 
        % Maximum conductance of I_NSCC
        g_NS = 0.0123; % [mS/uF]
        % Permeability of I_NSCC to Ca:Cs
        % defined as PnsCa=0.89 and PnsCs=1 separately in CellML; ratio is the same
        P_CaCs = 0.89; % [dimensionless]
        % Permeability of I_NSCC to Na:Cs
        % defined as PnsNa = 0.9 and PnsCs=1 separately in CellML; ratio is the same
        P_NaCs = 0.9; % [dimensionless]
        % Permeability of I_NSCC to K:Cs
        % defined as PnsK = 1.3 and PnsCs=1 separately in CellML; ratio is the same
        P_KCs = 1.3; % [dimensionless]
        % Half-saturation concentration for Mg_o inhibition of I_NSCC
        % Not explicitly defined in CellML, but used in I_NSCC equation as 0.281007
        K_dMg = 0.281007; % [mM] % Table S4: 0.28

        % Calcium-Activate Chloride Current (I_Cl(Ca))
        % Maximium conductance of I_Cl(Ca)
        g_Cl = 0.1875; % [mS/uF] - CellML

        % Na-K Exchanger Current (I_NaK)
        % Maximum conductance of I_NaK
        % ginak in CellML
        g_NaK = 1.7; % [A/F]
        % Half-saturation concentration for K_o dependency of I_NaK
        % nakKmko in CellML
        K_mK = 2.0; % [mM]
        % Hill coefficient for K_o dependency of I_NaK
        % not explicitly defined in CellML but used in knak equation as 1.5
        n_K = 1.5;
        % Half-saturation concentration for Na_i dependency of I_NaK
        % nakKmnai in CellML
        K_mNa = 22.0; % [mM]
        % Hill coefficient for Na_i dependency of I_NaK
        % not explicitly defined in CellML but used in nnak equation as 2
        n_Na = 2;

        % Calcium-Dependent Force
        % Half-saturation Ca concentration for force
        K_mF = 161.301 * 1e-6; % [mM]
        % Hill coefficient for Ca_i dependency of force
        n_F = 3.60205; % [dimensionless]
        % Maximal force
        F_max = 3.0 % [uN]

        % SR Ca2+ Channels/Fluxes - Shuai and Jung
        % v_IP3R = 0.6e-3; % [1/s] -> [1/ms] v_C 
        % % v_SERCA = 0.5e-3; % [1/s] -> [1/ms] v_P % should be mM/ms
        % % If mM/ms, then Src is correct
        % v_SERCA = 0.5e-6; % [uM/s] -> [mM/ms] v_P % should be mM/ms
        % 
        % k_SRleak = 0.001e-6; % [uM/s] -> [mM/ms] v_L
        % Ca_u = 15.0e-3; % [uM] -> [mM] Ca_ER
        % k_SERCA = 0.1e-3; % [uM] -> [mM] k
        % a = 0.2; % [1/(uM*s)] -> [1/(mM*ms)] (equal)
        % d_1 = 0.13e-3; % [uM] -> [mM]
        % d_2 = 1.05e-3; % [uM] -> [mM]
        % d_3 = 0.94e-3; % [uM] -> [mM]
        % d_4 = 0.13e-3; % [uM] -> [mM]
        % d_5 = 0.08e-3; % [uM] -> [mM]

        % SR Ca2+ Channels/Fluxes - Src Table
        % TODO: Double check units
        % v_IP3R = 1.95e-6; % [mM/ms] v_C
        % v_SERCA = 0.5e-6; % [mM/ms] v_P % CORRECT
        % k_SRleak = 0.001e-6; % [1/ms] v_L
        % Ca_u = 0.015; % [mM] Ca_ER
        % k_SERCA = 0.0001; % [mM] k
        % a = 0.2; % [1/(uM*s)] -> [1/(mM*ms)] (equal) not given
        % d_1 = 0.00105; % [mM/ms]
        % d_2 = 0.0005; % [mM/ms]
        % d_3 = 0.0003; % [mM/ms]
        % d_4 = 0.00013; % [mM/ms]
        % d_5 = 0.00008; % [mM]

        % It's possible d_1 and d_2 are swapped based on reading the src text - pg 12

        % CORRECT PARAMETERS IT SEEMS
        % SR Ca2+ Channels/Fluxes - src w/ Shuai and Jung d_i's + pg12 mods
        % TODO: Double check units
        v_IP3R = 1.95e-3; % [1/ms] v_C
        v_SERCA = 0.5e-6; % [mM/ms] v_P % CORRECT
        k_SRleak = 0.018e-3; % [1/ms] v_L
        Ca_u = 0.015; % [mM] Ca_ER
        k_SERCA = 0.0001; % [mM] k
        a = 0.2; % [1/(uM*s)] -> [1/(mM*ms)] (equal) not given
        d_1 = 0.0005; % [mM]
        d_2 = 0.00105; % [mM/ms]
        d_3 = 0.0003; % [mM/ms]
        d_4 = 0.00013; % [mM/ms]
        d_5 = 0.00008; % [mM]




        % % Rate of IP3 binding to IP3 receptors
        % v_IP3R = 1.95e-6; % [mM/ms] % v_C in Shuai and Jung
        % % Dissociation constant of the IP3-binding site on IP3R
        % % d_1 = 0.00105; % [mM/ms] - Src
        % % d_1 = 0.13e-3; % [mM] - Shuai & Jung
        % d_1 = 0.5e-3; % [mM] - Src Text, pg 12
        % % Scaling factor controlling
        % d_2 = 0.0005; % [mM/ms] - Src
        % % d_2 = 1.05e-3; % [mM] - Shuai & Jung
        % % IP3-related dissociation constant appearing
        % d_3 = 0.0003; % [mM/ms] - Src
        % % d_3 = 0.94e-3; % [mM] - Shuai & Jung
        % % Constant related to IP3 receptor activation
        % d_4 = 0.00013; % [mM/ms] - Src
        % % d_4 = 0.13e-3; % [mM] - Shuai & Jung
        % % Dissociation constant of the Ca2+ activation site on IP3R
        % d_5 = 0.00008; % [mM] - Src
        % % d_5 = 0.08e-3; % [mM] - Shuai & Jung
        % % Calcium concentration in the SR
        % Ca_u = 0.015; % [mM] - Src % good order of magnitude (comp to testrow's 0.023 mM from CellML)
        % % Ca_u = 15e-3; % [mM] - Shuai & Jung
        % % Maximum uptake rate of SERCA channel
        % v_SERCA = 0.5e-6; % [mM/ms] - src
        % % v_SERCA = 0.5e-3; % [1/s] - % v_P in Shuai & Jung
        % % The rate constant for calcium uptake via SERCA channel
        % k_SERCA = 0.0001; % [mM] % k in Shuai & Jung
        % % Rate constant for calcium leakage
        % k_SRleak = 1e-8; % [1/ms] % - src, v_L in Shuai & Jung
        % % k_SRleak = 0.001e-6; % [uM/s -> mM/ms] - Shuai & Jung
    end
end