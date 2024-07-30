function [Channel_est,Y_est] = esprit_dualwideband(Yn,MIMO_info,varargin)
% -------------------------------------------------------------------------
% esprit on each subband, save all the estimated AOD/AOAs


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
   fprintf('---------------esprit begins, on each subchannel---------------\n')
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

Channel_est.H = zeros(Nr,Nt,K);
Y_est = zeros(Q,P,K);

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
    
    A_ms_est = exp(   - 1i .* 2 .* pi .* (   1+K_select(k)*f_s/K_0/f_c   ) .* (0:Nr-1)' .* theta_est(k,:)   );
    A_bs_est = exp(   - 1i .* 2 .* pi .* (   1+K_select(k)*f_s/K_0/f_c   ) .* (0:Nt-1)' .* phi_est(k,:)   );
    W_times_Ams = W.' * A_ms_est;
    F_times_Abs = F.' * A_bs_est;
    Abs_khatri_Ams_est = khatrirao(F_times_Abs,W_times_Ams);
    alpha_aug_est(k,:) = (   (Abs_khatri_Ams_est'*Abs_khatri_Ams_est) \ (Abs_khatri_Ams_est'*Ynk(:))   ).';
    
    if par.show_info_during_process
        angleplot3(theta_est(k,:),phi_est(k,:),alpha_aug_est(k,:));
        hold on
    end
    
    Channel_est.H(:,:,k) = A_ms_est * diag(alpha_aug_est(k,:)) * A_bs_est.';
    Y_est(:,:,k) = W_times_Ams * diag(alpha_aug_est(k,:)) * F_times_Abs.';
end

Channel_est.theta = theta_est;
Channel_est.phi = phi_est;
Channel_est.alpha_aug = alpha_aug_est;

if par.show_info
    fprintf('esprit ends, time cost:%.2f,\n',toc)
end
end