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
% This function returns kron(B,A) from kron(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:  AkronB = kron(A,B)
%                               where  A: matrix of size NAh x NAv
%                                      B: matrix of size NBh x NBv
%
% Output Arguments: BkronA: kron(B,A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BkronA = switchKron(AkronB,NAh,NAv,NBh,NBv)
% horizontal indices rearrangement
h_ind = (1:NAh*NBh).';
H_ind = reshape(h_ind,[NBh,NAh]).';
h_ind_BA = H_ind(:);

% vertical indices rearrangement
v_ind = (1:NAv*NBv).';
V_ind = reshape(v_ind,[NBv,NAv]).';
v_ind_BA = V_ind(:);

BkronA = AkronB(h_ind_BA,v_ind_BA);

end

