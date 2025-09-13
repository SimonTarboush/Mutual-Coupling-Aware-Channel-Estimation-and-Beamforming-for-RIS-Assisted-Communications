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
% This function implements the OMP algorithm over separate multi-carriers (solve the compressed sensing problem for evey subcarrier)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% y: Measurements vector
% Upsilon: The sensing matrix 
% Pilot_pwr: The transmission power in the training phase
% Lbar: Estimation of the sparsity level
% K: Number of subcarriers
%
% Output Arguments:
% H_Est: The estimation of the channel in beamspace
% Support: The index selected by the OMP algorithm (known as the support)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H_Est, Support] = OMP_SMC(y,Upsilon,Pilot_pwr,Lbar,K)
R = y;
G = size(Upsilon,2);
H_Est = zeros(G,K);
Support = zeros(Lbar,K);
for indx_subc = 1:K
    % Solve the problem for every subcarrier
    Index = [];
    for ell =1: Lbar % threshold = sparsity level
        % Find AoD/AoA pair
        [~, indx_c] = max(abs((Upsilon'*R(:,indx_subc)).^2));
        % Update AoD/AoA pair index
        Index = [Index indx_c];
        h = zeros(G,1);
        % Estimate channel gains by solving the Least Square problem
        h(Index,:) = 1/sqrt(Pilot_pwr)*pinv(Upsilon(:,Index)'*Upsilon(:,Index))*Upsilon(:,Index)'*y(:,indx_subc);
        % Update residual
        R(:,indx_subc) = y(:,indx_subc) - sqrt(Pilot_pwr)*Upsilon(:,Index)*h(Index,:);
    end
    % Define the estimated channel support
    Support(:,indx_subc) = Index;
    % Define the estimated channel vector
    H_Est(:,indx_subc) = h;
end
end