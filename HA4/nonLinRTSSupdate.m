function [xs, Ps] = nonLinRTSSupdate(xs_kplus1, Ps_kplus1, xf_k, Pf_k, xp_kplus1, Pp_kplus1, ...
                                     f, T, sigmaPoints, type)

    n = length(xf_k);
    switch type
        case 'EKF'
            [~, Fx] = f(xf_k, T);
            Pk_kplus1 = Pf_k * Fx';
            
        case {'UKF', 'CKF'}

            [SP,W] = sigmaPoints(xf_k, Pf_k, type);
            
            Pk_kplus1 = zeros(n, n);
            for i = 1:numel(W)
                Pk_kplus1 = Pk_kplus1 + (SP(:,i) - xf_k)*( f(SP(:,i) ,T) - xp_kplus1)'*W(i);
            end
           
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end 

    G_k = Pk_kplus1/Pp_kplus1;
    xs = xf_k + G_k*(xs_kplus1 - xp_kplus1);
    Ps = Pf_k - G_k*(Pp_kplus1 - Ps_kplus1)*G_k';
end