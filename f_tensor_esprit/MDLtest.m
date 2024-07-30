function L_est = MDLtest(Y,MIMO_info,L_min)
% ------------------------------
% path number estimation by MDL
% ------------------------------
% written by XU Le
% ------------------------------
    P = MIMO_info.P;
    Q = MIMO_info.Q;
    K = MIMO_info.K;
    
    L_est = zeros(K,1);
    for kn = 1:K
        Y1 = Y(:,:,kn);
        e = eig(Y1*Y1'/P);
        e = flip(e); % descending order

        L_max = min(size(Y1));
        Vcomp = zeros(L_max,1);
        for i = 1:L_max
            numer = 1;
            domin = 0;
            for j = i+1:L_max
                numer = numer*e(j)^(1/(Q-i));
                domin = domin+e(j)/(Q-i);
            end
            if i ~= L_max
                Vcomp(i) = i/2*(2*Q-i)*log(P) - P*(Q-i)*log(numer/domin);
            else
                Vcomp(i) = i/2*(2*Q-i)*log(P);
            end

        end

        [~,L_est(kn)] = min(Vcomp);
        if L_est(kn)<L_min
            L_est(kn) = L_min;
        end
    end
    
    L_est = mode(L_est); % most frequent value
end