function [Channel_est,Y_est] = tensor_esprit_dualwideband(Yn,MIMO_info,varargin)
% -------------------------------------------------------------------------
% For the tensor esprit, the M2 strategy only composes a CPD tensor if
% there are multiple time slots (T), and P training pilots in each time
% slot. Furthermore, since W is fixed for each time slot, it is equivalent
% to open antenna 1-P for a time slot, then 1-P in the next round.

% Instead, we consider only one time slot, and the tensor esprit reduces to
% the traditional esprit for each subcarrier. As P-1>=L and Q>L, the
% classical ESPRIT should be unique ([]?).
% -------------------------------------------------------------------------
% Written by XU Le
% -------------------------------------------------------------------------

%%  parse parameters
p = inputParser;
addRequired(p,'Yn',@isnumeric);
addRequired(p,'MIMO_info',@isstruct);
addOptional(p,'is_path_number_known',false,@islogical);
addOptional(p,'default_path_number',4,@isnumeric);
addOptional(p,'show_info',true,@islogical);
addOptional(p,'show_info_during_process',false,@islogical);

parse(p,Yn,MIMO_info,varargin{:});
par = p.Results;

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

if par.show_info
   fprintf('---------------tensor esprit begins---------------\n')
   disp(par)
   tic
end

%% Estimate the number of path
if par.is_path_number_known
    L_est = par.default_path_number;
else
    min_path_number = 1;
    L_est = MDLtest(Yn,MIMO_info,min_path_number);
end

theta_est = zeros(K,L_est);
phi_est = zeros(K,L_est);
alpha_aug_est = zeros(K,L_est);

%% Channel parameter estimation for each subband using (tensor) esprit
for k = 1:K
    Ynk = Yn(:,:,k); % Q * P
    [U,Sigma,V] = svd(Ynk); % Ynk = U*Sigma*V'(hermitian)
    
    U = U(:,1:L_est); Sigma = Sigma(1:L_est,1:L_est); V = V(:,1:L_est);
    
    U1 = U(1:end-1,:); % (Q-1) * P
    U2 = U(2:end,:);
    
    MEMinv = (U1'*U1)\(U1'*U2); % P * P, with rank min(P,Q-1)
    [M,E] = eig(MEMinv);
    e = diag(E); % min(P,Q-1) * min(P,Q-1)
    
    % AOA estimation, as a column correspond to a received signal
    z_theta = e.' ./ abs(e.');
    theta_est(k,:) = - angle(z_theta)./2./pi./(   1+K_select(k)*f_s/K_0/f_c   );
    
    A_bs_subspace = conj(V) * Sigma.' * inv(M).'; % P * L_est
    
%     [U,Sigma,V] = svd(A_bs_subspace);
%     U = U(:,1:L_est); Sigma = Sigma(1:L_est,1:L_est); V = V(:,1:L_est);
%     MEMinv = (U1'*U1)\(U1'*U2); % P * P, with rank min(P,Q-1)
%     [M,E] = eig(MEMinv);
%     e = diag(E);
%     z_phi = e.' ./ abs(e.');
    
    
    z_phi = zeros(1,L_est);
    for ell = 1:L_est
        u1 = A_bs_subspace(1:end-1,ell);
        u2 = A_bs_subspace(2:end,ell);
        z_phi(ell) = u1'*u2 / (u1'*u1);
    end
    % AOD estimation
    phi_est(k,:) =  - angle(z_phi)./2./pi./(   1+K_select(k)*f_s/K_0/f_c   );
    [theta_est(k,:),index] = sort(theta_est(k,:));
    phi_est(k,:) = phi_est(k,index);
    
    if par.show_info_during_process
        angleplot3(theta_est(k,:),phi_est(k,:),zeros(L_est,1));
        hold on
    end
end

%% Refine the channel parameters based on all the subbands
% Channel_est.theta = mean(theta_est,1);
% Channel_est.phi = mean(phi_est,1);
[Channel_est.theta,Channel_est.phi] = anglepair(theta_est.',phi_est.',MIMO_info);
for k = 1:K
    ynk = reshape(   Yn(:,:,k), [Q*P,1]   ); % Q * P
    A_ms_est = exp(   - 1i .* 2 .* pi .* (   1+K_select(k)*f_s/K_0/f_c   ) .* (0:Nr-1)' .* Channel_est.theta   );
    A_bs_est = exp(   - 1i .* 2 .* pi .* (   1+K_select(k)*f_s/K_0/f_c   ) .* (0:Nt-1)' .* Channel_est.phi   );
    
    W_times_Ams = W.' * A_ms_est;
    F_times_Abs = F.' * A_bs_est;
    
    Abs_khatri_Ams_est = khatrirao(F_times_Abs,W_times_Ams);
    % alpha.*phaseshift estimation
    alpha_aug_est(k,:) = (   (Abs_khatri_Ams_est'*Abs_khatri_Ams_est) \ (Abs_khatri_Ams_est'*ynk)   ).';
end
Channel_est.tau = zeros(1,L_est);
Channel_est.alpha = zeros(1,L_est);
for ell = 1:L_est
    u1 = alpha_aug_est(1:end-1,ell);
    u2 = alpha_aug_est(2:end,ell);
    z_tau = u1'*u2 / (u1'*u1);
    Channel_est.tau(ell) = - angle(z_tau)/2/pi/(   ( K_select(2)-K_select(1) ) * f_s / K_0   );
    
    
%     tau_grid = 0:0.0001*1 / MIMO_info.f_s:1 / MIMO_info.f_s;
%     tau_dic = exp(   -1i * 2 * pi * K_select' * f_s / K_0 * tau_grid   );
%     
%     [~,max_index] = max(abs(tau_dic'*alpha_aug_est(:,ell)));
%     Channel_est.tau(ell) = tau_grid(max_index);
    
    g_tau = exp(   - 1i * 2 * pi * K_select'*f_s/K_0 * Channel_est.tau(ell)   ); % K * 1 phase shift vector in the K channels
    Channel_est.alpha(ell) = (g_tau'*g_tau) \ (g_tau'*alpha_aug_est(:,ell));
end

if par.show_info
    angleplot3(Channel_est.theta,Channel_est.phi,abs(Channel_est.alpha));
    fprintf('tensor esprit ends, time cost:%.2f,\n',toc)
end

%% Estimate the channel and the received signal
Channel_est.H = zeros(Nr,Nt,K); % only test the K channels, not K_0
Y_est = zeros(Q,P,K);
for k = 1:K
    A_ms_est = exp(   - 1i .* 2 .* pi .* (   1+K_select(k)*f_s/K_0/f_c   ) .* (0:Nr-1)' .* Channel_est.theta   );
    A_bs_est = exp(   - 1i .* 2 .* pi .* (   1+K_select(k)*f_s/K_0/f_c   ) .* (0:Nt-1)' .* Channel_est.phi   );
    g_tau = exp(   - 1i * 2 * pi * K_select(k)*f_s/K_0 * Channel_est.tau   );
    alpha_aug = Channel_est.alpha.*g_tau;
    
    Channel_est.H(:,:,k) = A_ms_est * diag(alpha_aug) * A_bs_est.';
    
    W_times_Ams = W.' * A_ms_est;
    F_times_Abs = F.' * A_bs_est;
    Y_est(:,:,k) = W_times_Ams * diag(alpha_aug) * F_times_Abs.';
end

end