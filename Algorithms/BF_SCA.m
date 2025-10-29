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
% This function implement the proposed MC-aware beamforming scheme
% following the successive convex approximation (SCA) presented in Sec. VI and
% Algorithm 2. The BS precoding and UE combining are presented in Sec. V-B
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
% Thres_SCA: the threshold on the change in the objective function between adjacent SCA iterations, below which the algorithm is considered to have converged
% init: indicator specifying whether to use random or predefined initialization
% combtype: indicator specifying whether to use analog or digital beamforming setup
%
% Output Arguments:
% w_opt: the optimized precoder at the BS
% f_opt: the optimized combiner at the UE
% gamma_opt: the optimized RIS response
% ga: the stored objective values across iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w_opt,f_opt,gamma_opt,ga] = BF_SCA(G,S,PI,PB,PU,NB,NU,Thres_Alt,Thres_SCA,init,combtype)

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
    
    [gamma,BFG_value] = SCA(S,q,B,gamma,PI,Thres_SCA,RPMC,G,w,f);
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


    %% Perform SCA procedure
    function [gamma_opt,BFG_value] = SCA(S,q,B,gamma,PI,Thres_SCA,RPMC,G,w,f)

        % Compute objective function and its derivatives
        objf = @(g,q,B) -g'*q*q'*g-g'*q*g.'*B*g-g'*B'*conj(g)*q'*g-g'*B'*conj(g)*g.'*B*g;
        df_g = @(g,q,B) -(g'*q)*conj(q)-(g'*q)*(B+B.')*g-(g'*B'*conj(g))*conj(q)-(g'*B'*conj(g))*(B+B.')*g;
        df_gc = @(g,q,B) -(q'*g)*q-(g.'*B*g)*q-(q'*g)*(B'+conj(B))*conj(g)-(g.'*B*g)*(B'+conj(B))*conj(g);

        N = size(S,1);
        go = gamma;
        ShrinkFactor = 0.1;
        ExpandFactor = 10;
        obj_old = +inf;
        obj_new = +inf;
        BFG_value = [];
        cnt1 = 0;

        while ~(abs(obj_new-obj_old) < Thres_SCA)

            % abs(obj_new-obj_old)

            % --- Construct surrogate function
            gi = go;
            gic = conj(go);
            Q = q*q' + q*gi.'*B + B'*gic*q' + B'*(gic*gi.')*B;
            [~,E] = eig(Q);
            eigv = sort(diag(E),'descend');
            K = eigv(1);
            b = ( K*gi' + gi'*q*gi.'*B.' + (gi'*B'*gic)*gi.'*B.' ).';
            c = ( K*gi.' + (q'*gi)*gi'*B' + (gi.'*B*gi)*gi'*B' ).';

            % --- Minimization surrogate function through KKT
            vlimit = [1,2];
            cnt2 = 0;
            while 1
                % Initial calculation
                Ps = abs((((K+vlimit(1))*eye(N)-Q.')^(-1)*b).'*(((K+vlimit(1))*eye(N)-Q)^(-1)*c));
                Pe = abs((((K+vlimit(2))*eye(N)-Q.')^(-1)*b).'*(((K+vlimit(2))*eye(N)-Q)^(-1)*c));
                % [Ps,Pe]
                % vlimit
                % search optimal v
                if abs(Ps-PI) < 1e-3 && abs(Pe-PI) < 1e-3
                    v = 0.5*vlimit*[1,1].';
                    break
                elseif abs(Ps-PI) < 5e-4
                    v = vlimit(1);
                    break
                elseif abs(Pe-PI) < 5e-4
                    v = vlimit(2);
                    break
                elseif Ps < PI && Pe < PI
                    vlimit(2) = vlimit(1);
                    vlimit(1) = vlimit(1)*ShrinkFactor;
                elseif Ps > PI && Pe > PI
                    vlimit(1) = vlimit(2);
                    vlimit(2) = vlimit(2)*ExpandFactor;
                elseif Ps > PI && Pe < PI
                    v = 0.5*vlimit*[1,1].';
                    Pn = norm((((K+v)*eye(N)-Q)^(-1)*c),2)^2;
                    switch Pn - PI < 0
                        case 1
                            vlimit(2) = v;
                        case 0
                            vlimit(1) = v;
                    end
                else
                    error('Unexpected errors!')
                end

                cnt2 = cnt2 + 1;
                if cnt2 > 100
                    break
                end
            end
            g = ((K+v)*eye(N)-Q)^(-1)*c;

            % --- Step-size line search to update the current feasible point
            alpha = 0.5;
            delta = 0.5;
            tk = 0;
            step = delta^tk;
            dir = g - go;
            f_go = objf(go,q,B);           % objective value at the current point
            grad_go = df_g(go,q,B);     % gradient w.r.t. \gamma
            grad_goc = df_gc(go,q,B);   % gradient w.r.t. conj(\gamma)
            while objf(go+step*dir,q,B) > f_go + alpha*step*(grad_go.'*dir+grad_goc.'*conj(dir))
                tk = tk + 1;
                step = delta^tk;
            end
            go = go + step*dir;

            % --- Evaluate cost function to decide if end the loop
            cost = -go'*( q*q' + q*go.'*B + B'*conj(go)*q' + B'*(conj(go)*go.')*B )*go;
            obj_old = obj_new;
            obj_new = cost;

            % --- save the beamforming gain at the current iteration
            BFG_current = cal_BFG(go,w,f,S,G,RPMC);
            BFG_value = [BFG_value,BFG_current];

            cnt1 = cnt1 + 1;
            if cnt1 > 5
                break
            end
        end

        gamma_opt = go;

    end

end

