%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Jointly developed by Simon Tarboush && Pinjun Zheng                          
% Last Modified: September, 2025
%
% If you use this code or any (modified) part of it in any publication, please cite the paper: 
% Pinjun Zheng*, Simon Tarboush*, Hadi Sarieddeen, and Tareq Y. Al-Naffouri, 
% "Mutual Coupling-Aware Channel Estimation and Beamforming for RIS-Assisted Communications", 
% accepted by IEEE Transactions on Wireless Communications.
% *: P. Zheng and S. Tarboush are co-first authors; they contributed equally to this paper and the implementation of these codes.
% IEEE-Xplore: https://ieeexplore.ieee.org/document/?????
% pre-print ArXiV: https://arxiv.org/abs/2410.04110
%
% Codes are also available on IEEE-Xplore: (will link the codes to CodeOcean)
% GitHub: https://github.com/SimonTarboush/Mutual-Coupling-Aware-Channel-Estimation-and-Beamforming-for-RIS-Assisted-Communications
%              https://github.com/ZPinjun/????
%
% Contact person email: simon.tarboush@tu-berlin.de && pinjun.zheng@ubc.ca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Key references:
% [1] Gradoni, Gabriele, and Marco Di Renzo. "End-to-end mutual coupling aware communication model for reconfigurable intelligent surfaces: An
% electromagnetic-compliant approach based on mutual impedances." IEEE Wireless Communications Letters 10.5 (2021): 938-942.   
% [2] R. Wang, Y. Yang, B. Makki and A. Shamim, "A Wideband Reconfigurable Intelligent Surface for 5G Millimeter-Wave Applications," in IEEE
% Transactions on Antennas and Propagation, vol. 72, no. 3, pp. 2399-2410, March 2024.   
% [3] P. Wang, J. Fang, H. Duan and H. Li, "Compressed Channel Estimation for Intelligent Reflecting Surface-Assisted Millimeter Wave Systems,"
% in IEEE Signal Processing Letters, vol. 27, pp. 905-909, 2020.   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates Fig. 8 in the paper 
% "Performance evaluation of spectral efficiency versus base station transmit power for different estimation and beamforming methods."
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Note that the results presented in the paper were obtained using larger
% dimensions for the RIS (N_I = 16×8), BS (N_B = 4×2), and UE (N_U = 2×1), 
% requiring substantial computational resources and a long simulation time. 
% Hence, in this script, smaller array dimensions are adopted to facilitate understanding and debugging the codes. 
%% Add Paths
path(pathdef); addpath(pwd);
cd Algorithms; addpath(genpath(pwd)); cd ..;
cd Utilities; addpath(genpath(pwd)); cd ..;
%%
clear;clc;
close all;
%% Evaluation setup
SaveSimulationResults = 1;              % Flag to indicate either (1) save the results and figures or (0) decline this operation
% parameters related to SNR
PU_prec_UL = 5;                         % Power constraint for uplink UE precoder (transmit power) in [mW] 
PB_comb_UL = 1;                         % Power constraint for uplink BS combiner in [mW] 
PB_DL_test = [-10,-5,0,5,10,15,20];     % Power constraint for downlink BS precoder (transmit power) in [dBm]
PU_DL = 1;                              % Power constraint for downlink UE combiner in [mW]
Amp_act_RIS = 7;                        % Average RIS amplification factors for RIS reflection vector power constraint
sigma2_tn = 10^(-95/10);                % Measured thermal noise power in [mW] (according to Fig. 25 of [2])
PB_DL_len = length(PB_DL_test);
% parameters related to number of trials
repeatNum = 10;                         % Number of Monte-Carlo simulation trials/runs/iterations
% parameters related to system dimensions (UE, BS, RIS)
arraySize_UE = [2,1];                   % UE array size, [horizontal dimension, vertical dimension] 
arraySize_BS = [3,1];                   % BS array size, [horizontal dimension, vertical dimension] 
arraySize_RIS = [8,6];                  % RIS array size, [horizontal dimension, vertical dimension]
% parameters related to mutual coupling at RIS
dr = 1/16;                              % RIS spacing [by lambda]
S_model = 'Analytical2';                % Selected method (options: 'Analytical1', 'Analytical2') to evaluate the mutual coupling expressed using S-parameters
% parameters related to proposed CS algorithm
M_B_rat = 0.75;                         % Ratio of the BS array size (less than 1) to express the number of measurement used in CS
M_I_rat = 0.75;                         % Ratio of the RIS array size (less than 1) to express the number of measurement used in CS
Lbar_BIU = 5;                           % Number of estimated cascaded paths by OMP (estimated sparsity since the exact number is unknown in practice)
DR_fact = 1;                            % Dictionary Reduction factor defined in Sec. IV-B3 and used in Algorithm 1 to controls how many atoms are kept after the DR algorithm 
% simulation methodology
TrainType = 'Variable';                 % Options: 'Static', 'Variable'. In 'Variable' (default), we use different training beams (at BS, UE, and RIS) for every trial
%% Set up System and Signals
% system and signal parameters
sp.fc = 30e9;                           % Center Frequency [Hz]
sp.c0 = 299792458;                      % Speed of Light [m/s]
sp.lambda_c = sp.c0/sp.fc;              % Wavelength
% Parameters for RIS, set the same as [1]    
sp.eta = 377;                           % The intrinsic impedance of free space [Ohm]
sp.k = 2*pi/sp.lambda_c;                % Wavenumber
sp.h = sp.lambda_c/64;                  % Antenna half-length 
sp.a = sp.lambda_c/500;                 % Antenna radius
sp.layout = 'YOZ';                      % The plane where the RIS is deployed on
% The system configuration, set the same as the far-field test setup in [2]     
% generate a UE
pos_UE = [-2.6*sind(30), 2.6*cosd(30), 0].';        % Position vector of UE in Cartesian coordinates
ori_UE = [-90, 0, 0].';                             % Euler angles, [0,0,0] corresponds to facing positive X-axis
spacing_UE = 0.5*sp.lambda_c;                       % UE element spacing
UE = c_array(pos_UE,ori_UE,arraySize_UE,spacing_UE);
% generate a BS
pos_BS = [2.2*sind(60), 2.2*cosd(60), 0].';       
ori_BS = [-90, 0, 0].';            
spacing_BS = 0.5*sp.lambda_c;   
BS = c_array(pos_BS,ori_BS,arraySize_BS,spacing_BS);
% generate a RIS
pos_RIS = [0, 0, 0].'; 
ori_RIS = [90, 0, 0].';    
spacing_RIS = dr*sp.lambda_c; 
RIS = c_array(pos_RIS,ori_RIS,arraySize_RIS,spacing_RIS);
% Visualize system layout
% UE.visualizeSys(BS,RIS); xlim([-3,3]); ylim([-1,5]);
% Compute scattering matrix for RIS mutual coupling 
S_mat = RIS.getS(sp,S_model); 
% Visualize S-parameters
% RIS.visualizeBP(10*log10(abs(S_mat))); title('RIS scattering matrix')
% Define system dimensions based on the generated instances of UE, BS, and RIS
N_U = prod(UE.dim);
N_B = prod(BS.dim);
N_I = prod(RIS.dim);
%% Far-field Uplink Channel Generation
% General parameters
L_U = 1;                        % Total number of paths between UE and RIS (LoS + NLoS)
PLE_IU = 2.1;                   % UE-RIS channel path loss exponent
L_B = 1;                        % Total number of paths between RIS and BS (LoS + NLoS)
PLE_BI = 2.2;                   % RIS-BS channel path loss exponent
% UE-RIS subchannel
chan_IU = c_channel(UE,RIS,L_U,PLE_IU,sp);
chan_IU.updateAoDAoA(0);
[H_IU,A_I_phi,A_U_varphi,Sigma_IU] = chan_IU.getH(0);    
% RIS-BS subchannel
chan_BI = c_channel(RIS,BS,L_B,PLE_BI,sp);
chan_BI.updateAoDAoA(0);
[H_BI,A_B_vartheta,A_I_theta,Sigma_BI] = chan_BI.getH(0);
% Ground truth equivalent cascaded channels
% Conventional cascaded UE-RIS-BS channel, check Sec. V-A and Eq. (31)
G_UL_cv = kron(A_U_varphi,A_B_vartheta)*kron(Sigma_IU,Sigma_BI)*get_KhatriRao_Prod(A_I_phi,A_I_theta).';
% MC-aware cascaded UE-RIS-BS channel, check Sec. V-A and Eq. (31)
G_UL_mc = kron(A_U_varphi,A_B_vartheta)*kron(Sigma_IU,Sigma_BI)*kron(A_I_phi,A_I_theta).';
%% Far-field Downlink Channel Generation
% BS-RIS subchannel
H_IB = H_BI.';
% RIS-UE subchannel
H_UI = H_IU.';
% Ground truth equivalent cascaded channels
% Conventional cascaded BS-RIS-UE channel, check Sec. V-A and Eq. (31)
G_DL_cv = kron(A_B_vartheta,A_U_varphi)*kron(Sigma_BI,Sigma_IU)*get_KhatriRao_Prod(A_I_theta,A_I_phi).';
% MC-aware cascaded BS-RIS-UE channel, check Sec. V-A and Eq. (31)
G_DL_mc = kron(A_B_vartheta,A_U_varphi)*kron(Sigma_BI,Sigma_IU)*kron(A_I_theta,A_I_phi).';
%% Construct the BS, UE, and RIS Dictionaries, Measurement, and Sensing matrices used for compressed-sensing (CS)-based channel estimation
G_ovs = 1;   % Grid oversampling (=1 : no oversampling in the angular domain)
G_I_h = ceil(G_ovs*RIS.dim(1));         % Quantization level of RIS horizontal dimension
G_I_v = ceil(G_ovs*RIS.dim(2));         % Quantization level of RIS vertical dimension
G_B_h = ceil(G_ovs*BS.dim(1));          % Quantization level of BS horizontal dimension
G_B_v = ceil(G_ovs*BS.dim(2));          % Quantization level of BS vertical dimension
G_U_h = ceil(G_ovs*UE.dim(1));          % Quantization level of UE horizontal dimension
G_U_v = ceil(G_ovs*UE.dim(2));          % Quantization level of UE vertical dimension
Abar_U_varphi = single(UE.getDictionary(G_U_h, G_U_v, UE.spacing/sp.lambda_c));
Abar_B_vartheta = single(BS.getDictionary(G_B_h, G_B_v, BS.spacing/sp.lambda_c));
Abar_UB = kron(Abar_U_varphi,Abar_B_vartheta);                  % following Eq. (24) in the paper
% RIS dictionary following conventional cascaded channel model 
Abar_I_phi = single(RIS.getDictionary(G_I_h, G_I_v, RIS.spacing/sp.lambda_c)); 
Abar_I_theta = single(RIS.getDictionary(G_I_h, G_I_v, RIS.spacing/sp.lambda_c));
Ahat_I_cv_pre = get_KhatriRao_Prod(Abar_I_phi,Abar_I_theta);    % following Eq. (24) in the paper
% Following [3, Proposition 1], the Khatri-Rao-based dictionary of conventional structure contains significant redundancy
Ahat_I_cv = unique(Ahat_I_cv_pre.','rows').';
% RIS dictionary following MC-aware cascaded channel model derived in the paper
Ahat_I = kron(Abar_I_phi,Abar_I_theta);                         % following Eq. (25) in the paper
% Full dictionary based on conventional cascaded channel model 
Dbar_cv = kron(Ahat_I_cv,Abar_UB);                              % following Eq. (24) in the paper
%% Precoder, Combiner, and RIS reflections Generation
seed_BF = 2;
P_I = (Amp_act_RIS^2)*N_I;                                      % Power constraint for RIS
% Traning beams at BS, UE, and RIS
M_B = ceil(N_B*M_B_rat);                                        % Number of training at BS
M_U = M_B;
M_I = ceil(N_I*M_I_rat);                                        % Number of training at RIS 
% Generate random precoding/combining weights at BS and UE
if strcmp(TrainType,'Static')
    W = BS.getRandBF(M_B,PB_comb_UL,seed_BF,'PhaseOnly');
    F = UE.getRandBF(M_U,PU_prec_UL,seed_BF,'PhaseOnly');
    P_norm = complex(zeros(M_B,N_B*N_U));
    for indx_mb = 1:M_B
        P_norm(indx_mb,:) = kron(transpose(F(:,indx_mb)),W(:,indx_mb)'); % following the definition in Eqs. (8) and (10) in the paper
    end
    % Generate random training weights (amplitude and phase) at RIS
    Theta_cv = RIS.getRandBF(M_I,P_I,seed_BF,'PhaseMag');       % following the definition in Eq. (16) in the paper
    Theta_mc = nan(N_I^2,M_I);
    for mi = 1:M_I
        Gamma = diag(Theta_cv(:,mi));
        Gamma_bar = (Gamma^(-1)-S_mat)^(-1);
        Theta_mc(:,mi) = Gamma_bar(:);                          % following the definition in Eq. (21) in the paper
    end
    % Measurement and Sensing matrices based on conventional cascaded channel model
    Psi_cv = kron(Theta_cv.',P_norm);                           % following the definition in Eq. (17) in the paper
    Xi_cv = Psi_cv*Dbar_cv;                                     % following the definition in Eq. (26) in the paper, case: cv
    % Measurement and Sensing matrices based on MC-aware cascaded channel model
    Psi_mc = kron(Theta_mc.',P_norm);                           % following the definition in Eq. (22) in the paper
end
%% Generate received signals (based on exact model, noise-free, with unit UE power)
if strcmp(TrainType,'Static')
    Y = zeros(M_B,M_I);
    for mi = 1:M_I
        Gamma = diag(Theta_cv(:,mi));
        Gamma_bar = (Gamma^(-1)-S_mat)^(-1);
        vecH_BIU = H_BI*Gamma_bar*H_IU;
        Y(:,mi) = P_norm*vecH_BIU(:);
    end
end
%% Evaluate spectral effciency based on estimated channels (fixed SNR) + beamforming (various SNRs)
getDLSNR = @(f,w,r,PU) abs(f'*H_UI*((diag(r))^(-1)-S_mat)^(-1)*H_IB*w)^2/(norm(f,2)^2*sigma2_tn + norm(f'*H_UI*((diag(r))^(-1)-S_mat)^(-1),2)^2*sigma2_tn);
% following the definition in Eq. (28) in the paper
getDLSE = @(f,w,r,PU) log2( 1 + getDLSNR(f,w,r,PU) ); % Check VII-B2

% In the paper we have presented results obtained by using digital combining (check Eq. (32) and the following discussion in Sec. V-B)
% Here, we provide the codes for both the analog and digital combining
% You should expect some differences with the obtained SE, however, the trend is similar
CombinerType = 'Digital';             % 'Analog', 'Digital'
data_GD_mc = nan(repeatNum,PB_DL_len);
data_SVD_mc = nan(repeatNum,PB_DL_len);
data_SCA_mc = nan(repeatNum,PB_DL_len);
data_GD_gt = nan(repeatNum,PB_DL_len);
data_SVD_gt = nan(repeatNum,PB_DL_len);
data_SCA_gt = nan(repeatNum,PB_DL_len);
tic;
for ind_rept = 1:repeatNum
    if strcmp(TrainType,'Static')
        % In static implementation once we change P_U that will change F and P = kron(F^T,W^H),
        % and since we do not need to generate a new realization from F and W,
        % we just adjust the power of P by the relevant factor
        P = sqrt(Psca_U)*P_norm;
        Y = zeros(M_B,M_I);
        for mi = 1:M_I
            Gamma = diag(Theta_cv(:,mi));
            Gamma_bar = (Gamma^(-1)-S_mat)^(-1);
            vecH_BIU = H_BI*Gamma_bar*H_IU;
            Y(:,mi) = P*vecH_BIU(:);
        end
    elseif strcmp(TrainType,'Variable')
        W = BS.getRandBF(M_B,PB_comb_UL,seed_BF,'PhaseOnly');
        F = UE.getRandBF(M_U,PU_prec_UL,seed_BF,'PhaseOnly');
        P = complex(zeros(M_B,N_B*N_U));
        for indx_mb = 1:M_B
            P(indx_mb,:) = kron(transpose(F(:,indx_mb)),W(:,indx_mb)');     % following the definition in Eqs. (8) and (10) in the paper
        end
        % Generate random training weights (amplitude and phase) at RIS
        Theta_cv = RIS.getRandBF(M_I,P_I,seed_BF,'PhaseMag');               % following the definition in Eq. (16) in the paper
        Theta_mc = nan(N_I^2,M_I);
        for mi = 1:M_I
            Gamma = diag(Theta_cv(:,mi));
            Gamma_bar = (Gamma^(-1)-S_mat)^(-1);
            Theta_mc(:,mi) = Gamma_bar(:);                                  % following the definition in Eq. (21) in the paper
        end
        Y = zeros(M_B,M_I);
        for mi = 1:M_I
            Gamma = diag(Theta_cv(:,mi));
            Gamma_bar = (Gamma^(-1)-S_mat)^(-1);
            vecH_BIU = H_BI*Gamma_bar*H_IU;
            Y(:,mi) = P*vecH_BIU(:);
        end
    else
        error('The only options are: ''Static'' and ''Variable''');
    end

    % Measurement and Sensing matrices based on conventional cascaded channel model
    Psi_cv = kron(Theta_cv.',P);    % following the definition in Eq. (17) in the paper
    Xi_cv = Psi_cv*Dbar_cv; % following the definition in Eq. (26) in the paper, case: cv
    % Measurement and Sensing matrices based on MC-aware cascaded channel model
    Psi_mc = kron(Theta_mc.',P);    % following the definition in Eq. (22) in the paper

    % Generate noise
    Noise = zeros(M_B,M_I);
    Noi_BS = randn(N_B,M_B,M_I)*sqrt(sigma2_tn/2) + 1j*randn(N_B,M_B,M_I)*sqrt(sigma2_tn/2);
    Noi_RIS = randn(N_I,M_B,M_I)*sqrt(sigma2_tn/2) + 1j*randn(N_I,M_B,M_I)*sqrt(sigma2_tn/2);
    for mi = 1:M_I
        Gamma = diag(Theta_cv(:,mi));
        Gamma_bar = (Gamma^(-1)-S_mat)^(-1);
        for mb = 1:M_B
            Noise(mb,mi) = W(:,mb)'*(Noi_BS(:,mb,mi) + H_BI*Gamma_bar*Noi_RIS(:,mb,mi)); % following Eq. (9) in the paper
        end
    end

    % Generate received signal
    Y_noi = Y + Noise;
    Y_GT = Y;

    % (1) Channel Estimation using conventional method (MC-unaware OMP)
    [sigma_UL_vec_cv, support_UL_cv] = OMP_SMC(Y_noi(:),Xi_cv,1,Lbar_BIU,1);
    Gest_UL_cv = reshape(Dbar_cv*sigma_UL_vec_cv,N_B*N_U,[]);

    % (2) Channel estimation using proposed Algorithm 1 in the paper "Two-Stage MC–aware with DR"
    % (2-a) Stage-1 "Coarse estimation" of the proposed algorithm is similar to (1) above, i.e., "MC-unaware OMP"
    % (2-b) Stage-2 "Refined estimation" use the prior information to construct a reduced size dictionary and then perform the estimation for MC-aware model

    % Get the initial estimate (prior information) from (1)
    sigma_UL_cv_mat =  reshape(sigma_UL_vec_cv,size(Abar_UB,2),size(Ahat_I_cv,2));
    [~, Supp_RIS_Est] = ind2sub(size(sigma_UL_cv_mat),support_UL_cv(1:Lbar_BIU));
    Abar_RISDict_Est = Ahat_I_cv(:,Supp_RIS_Est);                       % line 2 in Algorithm 1
    [~, indx_DR1] = get_Proposed_DR(Ahat_I(1:size(Ahat_I_cv,1),1:size(Ahat_I_cv,2)),Abar_RISDict_Est,size(Ahat_I_cv,2)*DR_fact(1)); % line 3 in Algorithm 1
    A_DR1 = Ahat_I(:,indx_DR1);
    Dbar_mc_DR1 = kron(A_DR1,Abar_UB);
    Xi_mc_DR1 = Psi_mc*Dbar_mc_DR1;                                     % line 4 in Algorithm 1
    % Perform DR-OMP for MC-aware estimation
    [sigma_UL_mc_DR1, ~] = OMP_SMC(Y_noi(:),Xi_mc_DR1,1,Lbar_BIU,1);    % line 5 in Algorithm 1
    Gest_UL_mc_DR = reshape(Dbar_mc_DR1*sigma_UL_mc_DR1,N_B*N_U,[]);

    % --- Conduct beamforming
    for ind_PBDL = 1:PB_DL_len
        disp(['RepeatNum: ',num2str(ind_rept),'/',num2str(repeatNum),'; P_B to be tested: [',num2str(PB_DL_test),']','; P_B: ', num2str(PB_DL_test(ind_PBDL))]);
        PB_DL = 10^(PB_DL_test(ind_PBDL)/10);

        % --- Beamforming using estimated channels
        Gest_DL_mc = switchKron(Gest_UL_mc_DR,N_U,N_I,N_B,N_I);
        % SCA optimization based on estimated mc channel
        [wSCA_est,fSCA_est,gammaSCA_est,gainsSCA_est] = BF_SCA(Gest_DL_mc,S_mat,P_I,PB_DL,PU_DL,N_B,N_U,1e-8,1e-2,'random',CombinerType);
        SE_SCA_est = getDLSE(fSCA_est,wSCA_est,gammaSCA_est,PB_DL);
        % GD optimization based on estimated mc channel
        [wGD_est,fGD_est,gammaGD_est,gainsGD_est] = BF_GD(Gest_DL_mc,S_mat,P_I,PB_DL,PU_DL,N_B,N_U,1e-8,1e-2,'random',CombinerType);
        SE_GD_est = getDLSE(fGD_est,wGD_est,gammaGD_est,PB_DL);
        % SVD optimization based on estimated cv channel
        Gest_DL_cv = switchKR(Gest_UL_cv,N_U,N_B);
        [wSVD_est,fSVD_est,gammaSVD_est,gainsSVD_est] = BF_SVD(Gest_DL_cv,S_mat,P_I,PB_DL,PU_DL,N_B,N_U,1e-8,1e-2,'random',CombinerType);
        SE_SVD_est = getDLSE(fSVD_est,wSVD_est,gammaSVD_est,PB_DL);

        % --- Beamforming using ground truth channels
        G_DL_mc = switchKron(G_UL_mc,N_U,N_I,N_B,N_I);
        % SCA optimization based on ground-truth mc channel
        [wSCA,fSCA,gammaSCA,gainsSCA_gt] = BF_SCA(G_DL_mc,S_mat,P_I,PB_DL,PU_DL,N_B,N_U,1e-8,1e-2,'random',CombinerType);
        SE_SCA = getDLSE(fSCA,wSCA,gammaSCA,PB_DL);
        % GD optimization based on ground-truth mc channel
        [wGD,fGD,gammaGD,gainsGD_gt] = BF_GD(G_DL_mc,S_mat,P_I,PB_DL,PU_DL,N_B,N_U,1e-8,1e-2,'random',CombinerType);
        SE_GD = getDLSE(fGD,wGD,gammaGD,PB_DL);
        % SVD optimization based on ground-truth cv channel
        G_DL_cv = switchKR(G_UL_cv,N_U,N_B);
        [wSVD,fSVD,gammaSVD,gainsSVD] = BF_SVD(G_DL_cv,S_mat,P_I,PB_DL,PU_DL,N_B,N_U,1e-8,1e-2,'random',CombinerType);
        SE_SVD = getDLSE(fSVD,wSVD,gammaSVD,PB_DL);

        data_GD_mc(ind_rept,ind_PBDL) = SE_GD_est;
        data_SVD_mc(ind_rept,ind_PBDL) = SE_SVD_est;
        data_SCA_mc(ind_rept,ind_PBDL) = SE_SCA_est;
        data_GD_gt(ind_rept,ind_PBDL) = SE_GD;
        data_SVD_gt(ind_rept,ind_PBDL) = SE_SVD;
        data_SCA_gt(ind_rept,ind_PBDL) = SE_SCA;
    end

end
Sim_Duration = toc;
hr = floor(Sim_Duration/3600);mint = floor((Sim_Duration - hr*3600)/60);sec = Sim_Duration - hr*3600 - mint*60;
fprintf('The simulation time is: %d hr %d min %f sec\n',hr,mint,sec);
%% Plot figure
figure
data = mean(data_SCA_gt);
plot(PB_DL_test,data,'--r'); hold on
data = mean(data_SCA_mc);
plot(PB_DL_test,data,'-rd'); hold on
data = mean(data_SVD_gt);
plot(PB_DL_test,data,'--b'); hold on
data = mean(data_SVD_mc);
plot(PB_DL_test,data,'-bs'); hold on
data = mean(data_GD_gt);
plot(PB_DL_test,data,'--g'); hold on
data = mean(data_GD_mc);
plot(PB_DL_test,data,'-go'); hold on
legend('Ground truth $\mathbf{G}_\mathrm{mc}$ + SCA BF', ...
    'Exact DR-OMP estimated $\hat{\mathbf{G}}_\mathrm{mc}$ + SCA BF', ...
    'Ground truth $\mathbf{G}_\mathrm{cv}$ + SVD BF', ...
    'Conv. improved OMP estimated $\hat{\mathbf{G}}_\mathrm{cv}$ + SVD BF', ...
    'Ground truth $\mathbf{G}_\mathrm{mc}$ + GD BF', ...
    'Exact DR-OMP estimated $\hat{\mathbf{G}}_\mathrm{mc}$ + GD BF', ...
    'interpreter','latex')   
xlabel('Transmit power (dBm)');
ylabel('Spectral Efficiency (bits/s/Hz)');
grid on
%%
if SaveSimulationResults
    cd Res
    simulationname = strcat(TrainType,'S',num2str(S_model),'SP',num2str(dr),'AMP',num2str(Amp_act_RIS),...
        'PUUL',num2str(PU_prec_UL),'PBUL',num2str(PB_comb_UL),'PUDL',num2str(PU_DL),...
        'I',num2str(repeatNum),'PBDLst',num2str(PB_DL_test(1)),'PBDLsp',num2str(PB_DL_test(end)),...
        'UR', num2str(L_U),'RB',num2str(L_B) , 'UE',num2str(arraySize_UE(1)),'x',num2str(arraySize_UE(2)),'BS',num2str(arraySize_BS(1)),'x',num2str(arraySize_BS(2)),...
        'RIS',num2str(arraySize_RIS(1)),'x',num2str(arraySize_RIS(2)),'M',num2str(M_B),'x',num2str(M_I),'OMP',num2str(Lbar_BIU),...
        'DR', num2str(DR_fact),num2str(int16(clock)));
    filename=strcat(simulationname,'.mat');
    save(sprintf('%s',filename),'PB_DL_test','sp','UE','BS','RIS','S_mat',...
        'data_GD_mc','data_SVD_mc','data_SCA_mc','data_GD_gt','data_SVD_gt','data_SCA_gt','-v7.3');
    GetfigHandles = findall(0,'Type','figure');
    for indx_fig = 1:numel(GetfigHandles)
        figname=strcat(simulationname,'_',num2str(indx_fig),'.fig');
        savefig(GetfigHandles(indx_fig),figname);
    end
    cd ..
end