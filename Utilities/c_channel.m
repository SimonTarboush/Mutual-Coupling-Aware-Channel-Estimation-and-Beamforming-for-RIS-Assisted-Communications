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
% This class generates the channel parameters and synthesizes the
% frequency-domain channel matrix for the multipath MIMO channel 
% between a pair of Tx and Rx given by two 'c_array' objects. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef c_channel < handle
    
    properties
        Tx       % Tx object 
        Rx      % Rx object
        L        % No. of paths (LoS + NLoS)
        PLE    % path loss exponent
        sp      % struct containing signal parameters
        AoD   % AoD at Tx, in [rad]
        AoA   % AoA at Rx, in [rad]
    end
    
    properties(Dependent)
        type    % channel type: 'SISO','SIMO','MISO','MIMO'
    end

    methods
        %% Constructor
        function obj = c_channel(Tx,Rx,L,PLE,sp,AoD,AoA)
            if nargin == 7
                obj.Tx = Tx;
                obj.Rx = Rx;
                obj.L = L;
                obj.PLE = PLE;
                obj.sp = sp;
                obj.AoD = AoD;
                obj.AoA = AoA;
            elseif nargin == 5
                obj.Tx = Tx;
                obj.Rx = Rx;
                obj.L = L;
                obj.PLE = PLE;
                obj.sp = sp;
            else
                error('Invalid inputs!')
            end
        end
        %% Get functions
        function type = get.type(obj)
            NT = prod(obj.Tx.dim);       % total number of antennas at Tx
            NR = prod(obj.Rx.dim);      % total number of antennas at Rx
            if NT == 1 && NR == 1          % SISO
                type = 'SISO';
            elseif NT == 1 && NR > 1      % SIMO
                type = 'SIMO';
            elseif NT > 1 && NR == 1      % MISO
                type = 'MISO';
            elseif NT > 1 && NR > 1        % MIMO
                type = 'MIMO';
            else
                error('Erroneous dimensions!')
            end
        end
        %% Regular Methods
        function updateAoDAoA(obj,seed) 
            % This method generates AoD and AoA between Tx and Rx based on deterministic LoS and random NLoS paths  
            rng(seed);
            switch obj.type
                case 'SISO'
                    obj.AoD = zeros(obj.L,2);
                    obj.AoA = zeros(obj.L,2);
                case 'SIMO'
                    Dist = norm(obj.Tx.pos-obj.Rx.pos,2);
                    obj.AoD = zeros(obj.L,2);
                    obj.AoA = zeros(obj.L,2);
                    % AoA
                    UnitDirVecLoc = obj.Rx.R.'*(obj.Tx.pos-obj.Rx.pos)/Dist;
                    Angle_azi = atan2(UnitDirVecLoc(2,:), UnitDirVecLoc(1,:));
                    Angle_ele = acos(UnitDirVecLoc(3,:));
                    obj.AoA(1,:) = [Angle_azi, Angle_ele];
                    obj.AoA(2:obj.L,1) = rand(obj.L-1,1)*pi-pi/2;
                    obj.AoA(2:obj.L,2) = rand(obj.L-1,1)*pi;
                case 'MISO'
                    Dist = norm(obj.Tx.pos-obj.Rx.pos,2);
                    obj.AoD = zeros(obj.L,2);
                    obj.AoA = zeros(obj.L,2);
                    % AoD
                    UnitDirVecLoc = obj.Tx.R.'*(obj.Rx.pos-obj.Tx.pos)/Dist;
                    Angle_azi = atan2(UnitDirVecLoc(2,:), UnitDirVecLoc(1,:));
                    Angle_ele = acos(UnitDirVecLoc(3,:));
                    obj.AoD(1,:) = [Angle_azi, Angle_ele];
                    obj.AoD(2:obj.L,1) = rand(obj.L-1,1)*pi-pi/2;
                    obj.AoD(2:obj.L,2) = rand(obj.L-1,1)*pi;
                case 'MIMO'
                    Dist = norm(obj.Tx.pos-obj.Rx.pos,2);
                    obj.AoD = zeros(obj.L,2);
                    obj.AoA = zeros(obj.L,2);
                    % AoD
                    UnitDirVecLoc = obj.Tx.R.'*(obj.Rx.pos-obj.Tx.pos)/Dist;
                    Angle_azi = atan2(UnitDirVecLoc(2,:), UnitDirVecLoc(1,:));
                    Angle_ele = acos(UnitDirVecLoc(3,:));
                    obj.AoD(1,:) = [Angle_azi, Angle_ele];
                    obj.AoD(2:obj.L,1) = rand(obj.L-1,1)*pi-pi/2;
                    obj.AoD(2:obj.L,2) = rand(obj.L-1,1)*pi;
                    % AoA
                    UnitDirVecLoc = obj.Rx.R.'*(obj.Tx.pos-obj.Rx.pos)/Dist;
                    Angle_azi = atan2(UnitDirVecLoc(2,:), UnitDirVecLoc(1,:));
                    Angle_ele = acos(UnitDirVecLoc(3,:));
                    obj.AoA(1,:) = [Angle_azi, Angle_ele];
                    obj.AoA(2:obj.L,1) = rand(obj.L-1,1)*pi-pi/2;
                    obj.AoA(2:obj.L,2) = rand(obj.L-1,1)*pi;
            end
        end
        function [H,A_R,A_T,Sigma] = getH(obj,seed)
            % This method generates channel matrix/vector/coeffcient between Tx and Rx 
            if exist('seed','var')
                rng(seed);
            end
            N = prod(obj.Tx.dim);
            M = prod(obj.Rx.dim);
            % compute spatial angles 
            SPang_Tx_h = (obj.Tx.spacing/obj.sp.lambda_c)*sin(obj.AoD(:,1)).*sin(obj.AoD(:,2));
            SPang_Tx_v = (obj.Tx.spacing/obj.sp.lambda_c)*cos(obj.AoD(:,2));
            SPang_Rx_h = (obj.Rx.spacing/obj.sp.lambda_c)*sin(obj.AoA(:,1)).*sin(obj.AoA(:,2));
            SPang_Rx_v = (obj.Rx.spacing/obj.sp.lambda_c)*cos(obj.AoA(:,2));
            % function to get array response
            A_oneaxis_resp = @(x,NumAnt) 1/sqrt(NumAnt)*exp(-1j*2*pi*(0:NumAnt-1)'*x);
            A_T = nan(N,obj.L);
            A_R = nan(M,obj.L);
            switch obj.type
                case 'SISO'
                    alpha = complex(zeros(obj.L,1));
                    alpha(1) = 1; % LoS
                    alpha(2:obj.L) = sqrt(1/2)*(randn(obj.L-1,1)+1j*randn(obj.L-1,1)); % L-1 NLoS
                    H = sum(alpha);
                    A_T = ones(1,obj.L);
                    A_R = ones(1,obj.L);
                case 'SIMO'
                    alpha = complex(zeros(obj.L,1));
                    alpha(1) = 1; % LoS
                    alpha(2:obj.L) = sqrt(1/2)*(randn(obj.L-1,1)+1j*randn(obj.L-1,1)); % L-1 NLoS
                    H = complex(zeros(M,N));
                    for indx_ell = 1:obj.L
                        a_R_y = A_oneaxis_resp(SPang_Rx_h(indx_ell),obj.Rx.dim(1));
                        a_R_z = A_oneaxis_resp(SPang_Rx_v(indx_ell),obj.Rx.dim(2));
                        a_R = kron(a_R_y,a_R_z);
                        H = H + alpha(indx_ell)*(a_R);
                        A_T(:,indx_ell) = 1;
                        A_R(:,indx_ell) = a_R;
                    end
                case 'MISO'
                    alpha = complex(zeros(obj.L,1));
                    alpha(1) = 1; % LoS
                    alpha(2:obj.L) = sqrt(1/2)*(randn(obj.L-1,1)+1j*randn(obj.L-1,1)); % L-1 NLoS
                    H = complex(zeros(M,N));
                    for indx_ell = 1:obj.L
                        a_T_y = A_oneaxis_resp(SPang_Tx_h(indx_ell),obj.Tx.dim(1));
                        a_T_z = A_oneaxis_resp(SPang_Tx_v(indx_ell),obj.Tx.dim(2));
                        a_T = kron(a_T_y,a_T_z);
                        H = H + alpha(indx_ell)*(a_T.');
                        A_T(:,indx_ell) = a_T;
                        A_R(:,indx_ell) = 1;
                    end
                case 'MIMO'
                    alpha = complex(zeros(obj.L,1));
                    alpha(1) = 1; % LoS
                    alpha(2:obj.L) = sqrt(1/2)*(randn(obj.L-1,1)+1j*randn(obj.L-1,1)); % L-1 NLoS
                    H = complex(zeros(M,N));
                    for indx_ell = 1:obj.L
                        a_T_y = A_oneaxis_resp(SPang_Tx_h(indx_ell),obj.Tx.dim(1));
                        a_T_z = A_oneaxis_resp(SPang_Tx_v(indx_ell),obj.Tx.dim(2));
                        a_T = kron(a_T_y,a_T_z);
                        a_R_y = A_oneaxis_resp(SPang_Rx_h(indx_ell),obj.Rx.dim(1));
                        a_R_z = A_oneaxis_resp(SPang_Rx_v(indx_ell),obj.Rx.dim(2));
                        a_R = kron(a_R_y,a_R_z);
                        H = H + alpha(indx_ell)*(a_R*a_T.');
                        A_T(:,indx_ell) = a_T;
                        A_R(:,indx_ell) = a_R; 
                    end
            end
            Dist = norm(obj.Tx.pos-obj.Rx.pos,2);
            H = (obj.sp.lambda_c/(4*pi*Dist))^(obj.PLE/2)*sqrt(N*M/obj.L)*H;
            Sigma = (obj.sp.lambda_c/(4*pi*Dist))^(obj.PLE/2)*sqrt(N*M/obj.L)*diag(alpha);
            % H_check = A_R*Sigma*A_T.';
        end
    end
end
