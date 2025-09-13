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
% This function computes the row-wise implementation of Khatri-Rao product of two matrices
%
% For two matrices with a similar row dimension the row-wise Khatri-Rao (KR) product
% is implemented by first transpose the input arrays, then perform the
% column-wise Khatri-Rao product and finally take the transpose of the
% output.
% the Column-wise KR is: Assume A of size (M1 x N) and B of size (M2 x N),
% the output is of size (M1M2 x N) and defined as
% KR_prod(A,B) = [kron(A(:,1),B(:,1)) .... kron(A(:,N),B(:,N))]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments: - A: matrix of size (M1 x N1)
%                              - B: matrix of size (M2 x N2)
% Output Arguments: C output matrix following row-wise Khatri-Rao product
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = get_KhatriRao_Prod(A, B)
% Transpose matrices for row-wise processing
A_T = A.';
B_T = B.';
[M_A_T, N_A_T] = size(A_T);
[M_B_T, N_B_T] = size(B_T);

if N_A_T ~= N_B_T
    error('Check dimensions: The two matrices Aand B should have the same number of rows !!')
end

C_T = zeros(M_A_T * M_B_T, N_A_T);
for indx_col=1:N_A_T
   % Perform Kronecker product on the columns of the transposed matrices
   col = kron(A_T(:,indx_col), B_T(:,indx_col));
   C_T(:,indx_col) = col(:);
end
% Transpose the result to get the row-wise Khatri-Rao product
C = C_T.';
end