%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Jointly developed by Simon Tarboush && Pinjun Zheng                          
% Last Modified: September, 2025
%
% If you use this code or any (modified) part of it in any publication, please cite the paper: 
% Pinjun Zheng*, Simon Tarboush*, Hadi Sarieddeen, and Tareq Y. Al-Naffouri, 
% "Mutual Coupling-Aware Channel Estimation and Beamforming for RIS-Assisted Communications", 
% accepted by IEEE Transactions on Wireless Communications.
% *: P. Zheng and S. Tarboush are co-first authors; they contributed equally to this paper and the implementation of these codes.
% IEEE-Xplore: https://ieeexplore.ieee.org/document/11176921
% pre-print ArXiV: https://arxiv.org/abs/2410.04110
%
% Codes are also available on IEEE-Xplore:
% GitHub: https://github.com/SimonTarboush/Mutual-Coupling-Aware-Channel-Estimation-and-Beamforming-for-RIS-Assisted-Communications
%              https://github.com/ZPinjun/MC-channel-estimation
%
% Contact person email: simon.tarboush@tu-berlin.de && pinjun.zheng@ubc.ca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function implements "DR-Algorithm" (lines 6, 7, and 8) in Algorithm 1 in the paper to construct the proposed Dictionary Reduction (DR) method
% The main aim is to select the most correlated columns from an (oversampled) dictionary.
% More details about Algorithm 1 in Sec. IV-B (DR-Algorithm: specifically Sec. IV-B3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% Abar_OVS: The (oversampled) dictionary
% Abar_Est: The selected columns of coarse estimation in Stage 1 of Algorithm 1
% G_DR: The desired size of the dictionary following the DR method
%
% Output Arguments:
% Abar_DR: The DR dictionary
% Indx_DR: The index of the columns selected from the "large" (oversampled) input dictionary with largest correlation value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Abar_DR, Indx_DR] = get_Proposed_DR(Abar_OVS, Abar_Est, G_DR)
G_DR = round(G_DR);
U_corr = Abar_OVS'*Abar_Est; % line 6 of Algorithm 1
[~ , sel_ind] = sort(sqrt(diag(U_corr*U_corr')),'descend'); 
Indx_DR = sel_ind(1:G_DR); % line 7 of Algorithm 1
Abar_DR = Abar_OVS(:,Indx_DR); % line 8 of Algorithm 1

end
