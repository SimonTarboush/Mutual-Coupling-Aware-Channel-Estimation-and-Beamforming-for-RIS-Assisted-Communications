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
% This function implement the baseline beamforming scheme
% following the SVD decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% G: the downlink exact equivalent cascaded channel
% S: the scattering matrix between RIS unit cells, unused in this function as SVD is a mutual coupling-unaware method 
% PI: the power of RIS response, i.e., constraint A in Eq. (30) in the paper
% PB: the transmit power at the BS
% PU: the combining power at the UE
% NB: the number of antennas at the BS
% NU: the number of antennas at the UE
% Thres_Alt: the threshold on the change in the objective function between adjacent alternating iterations, below which the algorithm is considered to have converged
% Thres_SCA: this variable is unused in this function 
% init: indicator specifying whether to use random or predefined initialization
% combtype: indicator specifying whether to use analog or digital beamforming setup
%
% Output Arguments:
% w_opt: the optimized precoder at the BS
% f_opt: the optimized combiner at the UE
% gamma_opt: the optimized RIS response
% ga: the stored objective values across iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w_opt,f_opt,gamma_opt,ga] = BF_SVD(G_cv,S,PI,PB,PU,NB,NU,Thres_Alt,Thres_SCA,init,combtype)

NI = size(G_cv,2);           % number of RIS elements

% Get initialization, random or designated by input
if strcmp(init,'random')
    rng(0);
    gamma = sqrt(PI/NI)*exp(1j*2*pi*rand(NI,1));
    w = sqrt(PB/NB)*exp(1j*2*pi*rand(NB,1));
    f = sqrt(PU/NU)*exp(1j*2*pi*rand(NU,1));
else
    gamma = init.gamma;
    w = init.w;
    f = init.f;
end

gain_old = cal_BFG(gamma,w,f,G_cv);
gain_new = inf;
ga.GainsAlt = [];
ga.GainsSCA = [];
ga.GainsAll = [];
count = 0;
gain_max = gain_old;
w_opt = w;
f_opt = f;
gamma_opt = gamma;

while ( (abs(gain_new-gain_old) > Thres_Alt) ) || ( (abs(gain_new-gain_old) < 1e-16) )

    if strcmp(combtype,'Analog')
        % Combining design (analog)
        H_UIB = reshape(G_cv*gamma,[NU,NB]);
        for ind = 1:10
            f = sqrt(PU/NU)*exp(1j*angle(H_UIB*w));
            w = sqrt(PB/NB)*exp(1j*angle(H_UIB'*f));
        end
        gain_old = gain_new;
        gain_new = cal_BFG(gamma,w,f,G_cv);
    else % Combining design (digital)
        H_UIB = reshape(G_cv*gamma,[NU,NB]);
        [U,~,V] = svd(H_UIB);
        w = sqrt(PB)*V(:,1);
        f = sqrt(PU)*U(:,1);
        gain_old = gain_new;
        gain_new = cal_BFG(gamma,w,f,G_cv);
    end
    % record the highest point in history
    if gain_new > gain_max
        w_opt = w;
        f_opt = f;
        gamma_opt = gamma;
        gain_max = gain_new;
    end

    % save objective function values
    ga.GainsAlt = [ga.GainsAlt,gain_new];
    ga.GainsSCA = [ga.GainsSCA,nan];
    ga.GainsAll = [ga.GainsAll,gain_new];

    % RIS configuration design
    t = G_cv.'*kron(w,conj(f));
    gamma = sqrt(PI/NI)*exp(-1j*angle(t)); 
    gain_old = gain_new;
    gain_new = cal_BFG(gamma,w,f,G_cv);

    % record the highest point in history
    if gain_new > gain_max
        w_opt = w;
        gamma_opt = gamma;
        gain_max = gain_new;
    end

    % save objective function values
    ga.GainsAlt = [ga.GainsAlt,nan];
    ga.GainsSCA = [ga.GainsSCA,gain_new];
    ga.GainsAll = [ga.GainsAll,gain_new];

    count = count + 1;
    if count > 20
        break
    end
end




    function gain = cal_BFG(gamma,w,f,G)
        gain = abs(kron(w.',f')*G*gamma(:));
        % disp(['Beamforming Gain: ', num2str(gain)]);
    end



end





