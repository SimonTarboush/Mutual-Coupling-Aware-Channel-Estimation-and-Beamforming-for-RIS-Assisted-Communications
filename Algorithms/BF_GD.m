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
% This function implement the MC-aware beamforming scheme following the
% gradient descent (GD) approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% G: the downlink exact equivalent cascaded channel
% S: the scattering matrix between RIS unit cells
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
function [w_opt,f_opt,gamma_opt,ga] = BF_GD(G,S,PI,PB,PU,NB,NU,Thres_Alt,Thres_SCA,init,combtype)

% Compute RIS response using accurate model
RPMC = @(gamma,S) (diag(gamma)^(-1) - S)^(-1);

NI = sqrt(size(G,2));           % number of RIS elements

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

gain_old = cal_BFG(gamma,w,f,S,G,RPMC);
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
        RP = RPMC(gamma,S);
        H_UIB = reshape(G*RP(:),[NU,NB]);
        for ind = 1:10
            f = sqrt(PU/NU)*exp(1j*angle(H_UIB*w));
            w = sqrt(PB/NB)*exp(1j*angle(H_UIB'*f));
        end
        gain_old = gain_new;
        gain_new = cal_BFG(gamma,w,f,S,G,RPMC);
    else % Combining design (digital)
        RP = RPMC(gamma,S);
        H_UIB = reshape(G*RP(:),[NU,NB]);
        [U,~,V] = svd(H_UIB);
        w = sqrt(PB)*V(:,1);
        f = sqrt(PU)*U(:,1);
        gain_old = gain_new;
        gain_new = cal_BFG(gamma,w,f,S,G,RPMC);
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
    t = G.'*kron(w,conj(f));
    q = diag(reshape(conj(t),[NI,NI]));
    B = S.*(reshape(t,[NI,NI])); 
    [gamma,BFG_value] = GD(gamma,q,B,S,PI,RPMC,G,w,f);
    gain_old = gain_new;
    gain_new = cal_BFG(gamma,w,f,S,G,RPMC);

    % record the highest point in history
    if gain_new > gain_max
        w_opt = w;
        gamma_opt = gamma;
        gain_max = gain_new;
    end

    % save objective function values
    ga.GainsAlt = [ga.GainsAlt,nan(size(BFG_value))];
    ga.GainsSCA = [ga.GainsSCA,BFG_value];
    ga.GainsAll = [ga.GainsAll,BFG_value];

    count = count + 1;
    if count > 20
        break
    end
end

    function gain = cal_BFG(gamma,w,f,S,G,RPMC)
        Gamma_bar = RPMC(gamma,S);
        gain = abs(kron(w.',f')*G*Gamma_bar(:));
        % disp(['Beamforming Gain: ', num2str(gain)]);
    end


    %% Perform GD procedure
    function [gamma_opt,BFG_value] = GD(gamma,q,B,S,PI,RPMC,G,w,f)

        % Compute objective function
        objf = @(g) -g'*q*q'*g-g'*q*g.'*B*g-g'*B'*conj(g)*q'*g-g'*B'*conj(g)*g.'*B*g;

        maxIter = 2;
        go = gamma;
        BFG_value = nan(1, maxIter);
        for iter = 1:maxIter

            % --- compute descent direction (according to Theorem 3.4 in [1])
            dir = ( q + (B+B.')*conj(go) ) * ( q'*go + go.'*B*go );
            dir = dir/norm(dir,2);
            
            % --- Step-size line search to update the current feasible
            % point (according to Algo. 9.2 in [2])
            alpha = 0.3;
            beta = 0.6;
            step = 1;
            dx = dir;
            % grad = -( q + (D+D.')*conj(go) ) * ( q'*go + go.'*D*go );
            grad = -( go'*q + go'*B*conj(go) ) * ( conj(q) + (B+B.')*go );
            cnt = 0;
            while 1
                f0 = objf(go);
                f1 = objf(go + step*dir);
                if f1 > f0 + alpha*step*grad'*dx
                    step = beta*step;
                else
                    break
                end
                cnt = cnt + 1;
                if cnt > 100
                    break
                end
            end
            go = go + step*dir;
            % --- project to maintain the power constraint
            go = sqrt(PI)*go/norm(go,2);

            % --- save the beamforming gain at the current iteration
            BFG_value(iter) = cal_BFG(go,w,f,S,G,RPMC);
        end
        gamma_opt = go;
    end

end





