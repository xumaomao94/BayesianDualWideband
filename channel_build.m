function [Yn,Ytrue,H] = channel_build(MIMO_info,Channel_info)

Nt = MIMO_info.Nt; % transmitting antennas
P = MIMO_info.P;
Nr = MIMO_info.Nr; % receiving antennas
Q = MIMO_info.Q;

K_0 = MIMO_info.K_0; % total number of subcarriers
K = MIMO_info.K; % number of used subcarriers
K_select = MIMO_info.K_select;

f_c = MIMO_info.f_c; % carrier frequency
f_s = MIMO_info.f_s; % subcarrier frequency

F = MIMO_info.F; % Nt*P
W = MIMO_info.W; % Nr*Q


L = Channel_info.L;  % number of path
SNR = Channel_info.SNR; % \dB
alpha = Channel_info.alpha; %size: L*1
tau = Channel_info.tau;
phi_eq = sin(Channel_info.phi)./2; % antenna space: lambda_c/2
theta_eq = sin(Channel_info.theta)./2;

% steering matrix
A_bs = zeros(Nt,L,K);
A_ms = zeros(Nr,L,K);
H = zeros(Nr,Nt,K);
Ytrue = zeros(Q,P,K);
Yn = zeros(Q,P,K);


for k = 1:K
    f_k = K_select(k) * f_s / K_0;
    A_bs(:,:,k) = exp(   -1i.*2.*pi .* (1+f_k/f_c) .* (0:Nt-1)' .* phi_eq'   );
    A_ms(:,:,k) = exp(   -1i.*2.*pi .* (1+f_k/f_c) .* (0:Nr-1)' .* theta_eq'   );
    
    alpha_aug = alpha .* exp(   -1i.*2.*pi .* f_k .* tau   );
    H(:,:,k) = A_ms(:,:,k) * diag(alpha_aug) * A_bs(:,:,k).'; % Nr * Nt
    if size(F,3) == 1
        Yk = H(:,:,k) * F; % Nr * P
    elseif size(F,3) == K
        Yk = H(:,:,k) * F(:,:,k);
    end
    
    Npower = 1 / 10^(SNR/10);
    W_rest = W(Q+1:end,:);
    if sum(abs(W_rest(:))) == 0 % Identity
        Nk = sqrt(Npower/Nr) .* sqrt(1/2) .* (   randn(Nr,P) + 1i*randn(Nr,P)   );
    else
        Nk = sqrt(Npower) .* sqrt(1/2) .* (   randn(Nr,P) + 1i*randn(Nr,P)   );
    end

    Ykn = Yk + Nk;
    if size(W,3) == 1
        Ytrue(:,:,k) = W.' * Yk;
        Yn(:,:,k) = W.' * Ykn;
    elseif size(W,3) == K
        Ytrue(:,:,k) = W(:,:,k).' * Yk;
        Yn(:,:,k) = W(:,:,k).' * Ykn;
    end
end

end