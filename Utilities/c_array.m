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
% This class generates an object that might represent a UE, BS, or RIS
% array at a specific location (based on a position and orientation). It
% also generates the scattering parameters and some essential functions
% related to create far-field (oversampled) dictionary and random beamforming weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef c_array < handle
    
    properties
        pos = [0, 0, 0].';             % array position [meter]
        ori = [30, 60, 90].';        % array orientation [deg]
        dim = [1, 1];                  % array size, [h dim, v dim]
        spacing = 0.005;           % array element spacing [meter]
        Gamma                         % array configuration
    end

    properties(Dependent)
        R       % rotation matrix of RIS
        P_i     % coordinates of RIS elements in body coordinate system (3xN) 
    end
    
    methods
        %% Constructor
        function obj = c_array(pos,ori,size,spacing)
            if nargin > 1           % antenna array
                obj.pos = pos;
                obj.ori = ori;
                obj.dim = size;
                obj.spacing = spacing;
            elseif nargin == 1  % single antenna
                obj.pos = pos;
            end
        end
        %% Get functions
        function R = get.R(obj)
            R = eul2rotm(reshape(deg2rad(obj.ori), [1,3]), 'ZYX');
        end
        function P_i = get.P_i(obj)
            yrange = (  (0:obj.dim(1)-1) - 0.5*(obj.dim(1)-1)  ) * obj.spacing;
            zrange = (  (0:obj.dim(2)-1) - 0.5*(obj.dim(2)-1)  ) * obj.spacing;
            P_i = [zeros(1,obj.dim(1)*obj.dim(2));
                kron(yrange, ones(1,obj.dim(2)));
                kron(ones(1,obj.dim(1)), zrange)];
        end
        %% Regular Methods
        function visualizeBP(~,BP)
            % Visualize a input pattern
            figure();imagesc(BP);colorbar;
        end
        function visualizeSys(obj,Rx,Reflector)
            % Visualiza the system layout. Obj is set as Tx.   
            UE = obj;
            BS = Rx;
            RIS = Reflector;
            figure
            % plot UE
            possr = get_array_shape(UE,0.4);
            possr = possr + UE.pos(1:2)*ones(1,size(possr,2));
            scatter(possr(1,:),possr(2,:),'xk','LineWidth',1.5); hold on
            % plot BS
            possr = get_array_shape(BS,0.4);
            possr = possr + BS.pos(1:2)*ones(1,size(possr,2));
            scatter(possr(1,:),possr(2,:),'b+','LineWidth',1.5); hold on
            % plot RIS
            possr = get_array_shape(RIS,0.3);
            possr = possr + RIS.pos(1:2)*ones(1,size(possr,2));
            scatter(possr(1,:),possr(2,:),'r.','LineWidth',1.5); hold on
            % plot path
            quiver(UE.pos(1),UE.pos(2),RIS.pos(1)-UE.pos(1),RIS.pos(2)-UE.pos(2),0); hold on
            quiver(RIS.pos(1),RIS.pos(2),BS.pos(1)-RIS.pos(1),BS.pos(2)-RIS.pos(2),0); hold on
            xlim([-50,50]);
            ylim([-50,50]);
            title('System layout');
            xlabel('x (m)');
            ylabel('y (m)');
            function possr = get_array_shape(array,scale)
                orientation = array.ori;
                length = array.dim(1);
                % a RIS pattern
                poss = [zeros(1,length);linspace(-scale*length/2,scale*length/2,length)];
                pos0 = [0.2;0];
                % rotate RIS patten
                gamma = orientation(1);
                possr = [-sind(gamma)*poss(2,:); cosd(gamma)*poss(2,:)];
                pos0r = [cosd(gamma)*pos0(1); sind(gamma)*pos0(1)];
                possr = [possr,pos0r];
            end
        end
        function S = getS(obj, sp, model)
            % Inputs:
            % - P: a 3xN matrix storing the positions of antennas, P = [p1,p2,...,pN].
            % - sp: a struct containing necessary antenna parameters
            % - model: choose which model to generate S
            % Output:
            % - S: a NxN matrix containing S parameters between each pair of antennas

            % ------ 
            % References:
            % [1] G. Gradoni, and M. Di Renzo. "End-to-end mutual coupling aware communication model for reconfigurable intelligent surfaces: An
            % electromagnetic-compliant approach based on mutual impedances," in IEEE Wireless Communications Letters 10.5 (2021): 938-942.

            % [2] P. Zheng, R. Wang, A. Shamim, and T. Y. Al-Naffouri, "Mutual Coupling in RIS-Aided Communication: Model Training and Experimental Validation," 
            % in IEEE Transactions on Wireless Communications, vol. 23, no. 11, pp. 17174-17188, Nov. 2024.

            % [3] A. Abrardo, A. Toccafondi, and M. Di Renzo. "Design of Reconfigurable Intelligent Surfaces by Using S-Parameter Multiport Network Theory—Optimization and Full-Wave Validation," 
            % in IEEE Transactions on Wireless Communications, vol. 23, no. 11, pp. 17084-17102, Nov. 2024

            % Adjust coordinates based on the antenna array layout. The path integrals
            % for the mutual/self impedance calculation in 'func_MutuImp_antenna.m' are
            % performed along the z-axis. So we need to shift coordinates to make sure
            % z-axis is the normal direction of the antenna array.
            P = obj.P_i;
            switch sp.layout
                case 'XOY'
                    % not implemented
                case 'XOZ'
                    P = [P(1,:);P(3,:);P(2,:)];
                case 'YOZ'
                    P = [P(2,:);P(3,:);P(1,:)];
                otherwise
                    disp('Unexcepted array layout!');
            end

            N = size(P,2);
            Z0 = 50;    % the reference impedance, typically 50 [Ohm]

            % Generate scattering matrix S based on different models
            % - 'Analytical1': Calculate mutual & self impedances (Z-parameters)
            %   according to [1], and then transform to S-parameters. This method
            %   usually overestimates the self impedances.
            % - 'Analytical2': Calculate mutual impedances (Z-parameters) according to
            %   [1], while assign self impedances as reference impedance (50 Ohm),
            %   i.e., assume no self impedances. Then, transform to S-parameetrs.
            % - 'Measured': Assign S-parameters based on the real measurement in [2].
            switch model
                case 'Analytical1'
                    Z = zeros(N,N);
                    % Obtain self impedances
                    SelfImp = get_MutuImp_antenna(P(:,1), P(:,1), sp);
                    % Obtain mutual impedances
                    for i = 1:N
                        p_i = P(:,i);
                        for j = i+1:N
                            p_j = P(:,j);
                            Z(i,j) = get_MutuImp_antenna(p_i, p_j, sp);
                            Z(j,i) = Z(i,j);
                        end
                    end
                    Z = Z +SelfImp*eye(N);
                    % Transform to S-parameters, based on [3, Eq.(24)]
                    S = (Z + Z0*eye(N))^(-1)*(Z - Z0*eye(N));
                case 'Analytical2'
                    Z = zeros(N,N);
                    % Set self impedances equal to Z0, based on [3, Eq.(26)]
                    SelfImp = Z0;
                    % Obtain mutual impedances
                    for i = 1:N
                        p_i = P(:,i);
                        for j = i+1:N
                            p_j = P(:,j);
                            Z(i,j) = get_MutuImp_antenna(p_i, p_j, sp);
                            Z(j,i) = Z(i,j);
                        end
                    end
                    Z = Z +SelfImp*eye(N);
                    % Transform to S-parameters, based on [3, Eq.(24)]
                    S = (Z + Z0*eye(N))^(-1)*(Z - Z0*eye(N));
                otherwise
                    warning('Unexpected model type!')
            end

            % This function computes the mutual impedance of two antennas according to Eq. (2) & (3) of [2]
            % ------
            % References:
            % [1] Gradoni, Gabriele, and Marco Di Renzo. "End-to-end mutual coupling aware communication model for reconfigurable intelligent surfaces: An electromagnetic-compliant approach based on mutual impedances." IEEE Wireless Communications Letters 10.5 (2021): 938-942.
            % [2] Di Renzo, Marco, Vincenzo Galdi, and Giuseppe Castaldi. "Modeling the Mutual Coupling of Reconfigurable Metasurfaces." 17th European Conference on Antennas and Propagation (EuCAP). IEEE, 2023.
            % ------
            function z_qp = get_MutuImp_antenna(pp, pq, sp)
                % ----- Input:
                % pp: center position of the antenna p
                % pq: center position of the antenna q
                % sp: other parameters and constants
                % ----- Output:
                % z_qp: mutual/self impedance between antennas p & q

                % get parameters
                eta = sp.eta;       % the intrinsic impedance of free space
                k0 = sp.k;          % wavenumber
                hp = sp.h;          % half-length of antenna p
                hq = sp.h;          % half-length of antenna q
                aq = sp.a;          % antenna radius

                % some preparatory calculations
                if norm(pp-pq,2) == 0       % if p = q
                    rho1 = aq;
                    rho2 = 0;
                else                        % if p != q
                    rho1 = sqrt( (pp(1) - pq(1))^2 + (pp(2) - pq(2))^2 );
                    rho2 = pp(3) - pq(3);
                end

                % compute integral Eq. (2) of [2]
                R1  = @(xi,z) sqrt(rho1^2 + (z - xi + rho2).^2);
                f1 = @(xi,z) exp(-1j*k0*R1(xi,z)).*sin(k0*(hp-abs(xi))).*sin(k0*(hq-abs(z))) ./ ( R1(xi,z)*sin(k0*hp)*sin(k0*hq) );
                f2 = @(xi,z) k0^2 - 1j*k0./R1(xi,z) - (k0^2*(z-xi+rho2).^2+1)./R1(xi,z).^2 + 3*1j*k0*(z-xi+rho2).^2./R1(xi,z).^3 + 3*(z-xi+rho2).^2./R1(xi,z).^4;
                f  = @(xi,z) f1(xi,z) .* f2(xi,z);
                z_qp = 1j*eta/(4*pi*k0) * integral2(f, -hp, hp, -hq, hq);
            end
        end
        function A_UPA = getDictionary(obj, G_h, G_v, d_ae)
            % Inputs:
            % obj:      Array object
            % G_h:      Grid size over horizontal dimension to quantize the angular region
            % G_v:      Grid size over vertical dimension to quantize the angular region
            % d_ae:     normalized antenna elements spacing
            % Output:
            % A:        The desired dictionary

            % Spatial Angles
            horz_spang = -1+1/G_h:2/G_h:1-1/G_h;
            vert_spang = -1+1/G_v:2/G_v:1-1/G_v;
            % compute array response vectors
            A_ULA_resp = @(x,NumAnt) 1/sqrt(NumAnt)*exp(-1j*2*pi*(0:NumAnt-1)'*x);
            A_h = A_ULA_resp(d_ae*horz_spang,obj.dim(1));  % A_h discretize the azimuth plane
            A_v = A_ULA_resp(d_ae*vert_spang,obj.dim(2));   % A_v discretize the elevation plane
            A_UPA = kron(A_h,A_v);
        end
        function W = getRandBF(obj,M,P,seed,option)
            % Inputs:
            % obj:           Array object
            % M:             Numer of beamformers
            % P:              Power constraint for each beamformer
            % seed:        Random number generator seed value
            % option:     'PhaseOnly': this indicates that we only change the phase shifter of the precoder/combiner/ ... and assume constant magnitude modulus. 
            %                  'PhaseMag': the second option is to change both the amplitude and phase of the precoder/combiner/ ... 
            %                  
            % Output:
            % W:        Random beamformer codebook, containing M columns.
            % Each column of W maintains the same power ||W(:,i)||^2 = P.
            
            if exist('seed','var')
                rng(seed);
            end
            N = prod(obj.dim);
            phase_rand = 2*pi*rand(N,M);
            if strcmp(option,'PhaseOnly')
                magnitude = ones(N,M);
            elseif strcmp(option,'PhaseMag')
                magnitude = rand(N,M);
            else
                error('The only options are: ''PhaseOnly'' and ''PhaseMag''');
            end
            scale = ones(N,1)*sqrt(ones(1,N)*(magnitude.*magnitude));
            magnitude = (magnitude./scale)*sqrt(P);
            W = magnitude.*exp(1j*phase_rand);
        end
    end
end