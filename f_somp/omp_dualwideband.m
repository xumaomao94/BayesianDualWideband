function [Channel_est,Y_est] = omp_dualwideband(Yn,MIMO_info,varargin)
% -------------------------------------------------------------------------
% Since the dictionaries are different under different channel frequencies,
% the omp algorithm is adopted to estimate each subchannel.

% Then Channel_est.H \in C^{Nr*Nt*K} is built for each k = 1 to K. And the
% AOD/AOA/alpha are different for different frequencies.
% -------------------------------------------------------------------------
% Written by XU Le
% -------------------------------------------------------------------------

%%  parse parameters
p = inputParser;
addRequired(p,'Yn',@isnumeric);
addRequired(p,'MIMO_info',@isstruct);
addOptional(p,'grid_sample_method',"sinspace",@isstring); % anglespace: -pi/2 ~ pi/2; sinspace: -0.5 ~ 0.5
addOptional(p,'dense_phi',128,@isnumeric);
addOptional(p,'dense_theta',128,@isnumeric);
addOptional(p,'max_path_number',4,@isnumeric);
addOptional(p,'tol',1e-2,@isnumeric);
addOptional(p,'whitening_or_not',false,@islogical);
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

tol = par.tol;
dense_phi = par.dense_phi;
dense_theta = par.dense_theta;
dense_grid = dense_phi*dense_theta;

if par.show_info
   fprintf('---------------somp begins---------------\n')
   disp(par)
   tic
end

% whitening
if par.whitening_or_not
    [V,D] = eig(MIMO_info.W.'*conj(MIMO_info.W));
    whiteningM = diag(1./sqrt(diag( D )))*V';
    Yn = whiteningM * reshape(Yn,MIMO_info.Q,[]);
    Yn = reshape(Yn,[MIMO_info.Q,MIMO_info.P,MIMO_info.K]);
    MIMO_info.W = MIMO_info.W * whiteningM.';
end


%% Channel estimation through omp
Channel_est.H = zeros(Nr,Nt,K);
Channel_est.path_exist_or_not = zeros(dense_grid,K);
Channel_est.phi = zeros(dense_grid,K);
Channel_est.theta = zeros(dense_grid,K);
Channel_est.alpha_aug = zeros(dense_grid,K);
Y_est = zeros(Q,P,K);

for k = 1:K
    [Dic_k,Dic_sqr_k] = dic_init_omp(MIMO_info,dense_phi,dense_theta,k,par.grid_sample_method);
    ynk = reshape(   Yn(:,:,k),[Q*P,1]   );
    res = ynk;
    
    Basis = [];
    ind_select = [];
    alpha_aug = [];
    
    % omp for the k-th channel
    for r = 1:min(   max(size(Dic_k,2)),par.max_path_number   )
        weight = abs(   Dic_k' * res    )./ Dic_sqr_k;%./ Dic_sqr_k;
        
        [~,ind] = max(weight);
        ind = ind(1); % in case of two max values
        x = Dic_k(:,ind)'*res / Dic_sqr_k(ind);%/ (   Dic_k(:,ind)'*Dic_k(:,ind)   );%
        res = res - Dic_k(:,ind)*x;
        
        if abs(x) < tol
            if r > 1
                break
            end
        end
        
        ind_select = [ind_select,ind];
        Basis = [Basis, Dic_k(:,ind)];
        Dic_k(:,ind) = 0;
        Dic_sqr_k(ind) = 1e20;
        alpha_aug = [alpha_aug,x];
        
        if par.show_info_during_process
            iphi = ceil(ind/dense_theta);
            itheta = ind - (iphi-1).*dense_theta;

            if strcmp(par.grid_sample_method,"anglespace")
                phi_select = -pi/2 + pi/2/dense_phi + pi/dense_phi.*(iphi-1);
                theta_select = -pi/2 + pi/2/dense_theta + pi/dense_theta.*(itheta-1);
                angleplot3(   sin(theta_select)/2,sin(phi_select)/2,abs(x)   );
            elseif strcmp(par.grid_sample_method,"sinspace")
                phi_eq_select = -0.5 + 0.5/dense_phi + 1/dense_phi.*(iphi-1);
                theta_eq_select = -0.5 + 0.5/dense_theta + 1/dense_theta.*(itheta-1);
                angleplot3(   theta_eq_select,phi_eq_select,abs(x)   );
            end
            hold on;
        end
    end
    
    % the corresponding equivalent angles
    Channel_est.path_exist_or_not(ind_select,k) = 1;
    iphi = ceil(ind_select/dense_theta);
    itheta = ind_select - (iphi-1).*dense_theta;
    if strcmp(par.grid_sample_method,"anglespace")
        phi_select = -pi/2 + pi/2/dense_phi + pi/dense_phi.*(iphi-1);
        theta_select = -pi/2 + pi/2/dense_theta + pi/dense_theta.*(itheta-1);
        Channel_est.phi(ind_select,k) = sin(phi_select)./2;
        Channel_est.theta(ind_select,k) = sin(theta_select)./2;
    elseif strcmp(par.grid_sample_method,"sinspace")
        Channel_est.phi(ind_select,k) = -0.5 + 0.5/dense_phi + 1/dense_phi.*(iphi-1);
        Channel_est.theta(ind_select,k) = -0.5 + 0.5/dense_theta + 1/dense_theta.*(itheta-1);
    end
    if par.show_info
%         if k == 1 % draw only first of the subbands for now
            angleplot3(   Channel_est.theta(ind_select,1),Channel_est.phi(ind_select,1),abs(alpha_aug)   );
            hold on;
%         end
    end
    
    % rebuild the channel and received signal
    f_k = K_select(k)*f_s/K_0;
    alpha_aug = (Basis'*Basis) \ (Basis'*ynk);
    Channel_est.alpha_aug(ind_select,k) = alpha_aug;
    Channel_est.H(:,:,k) = exp(   -1i.*2.*pi .* (1+f_k/f_c) .* (0:Nr-1)' * Channel_est.theta(ind_select,k).'   )...
        * diag(alpha_aug) * exp(   -1i.*2.*pi .* (1+f_k/f_c) .* (0:Nt-1)' * Channel_est.phi(ind_select,k).'   ).';
    Y_est(:,:,k) = MIMO_info.W.' * Channel_est.H(:,:,k) * F;
    if par.whitening_or_not
        Y_est(:,:,k) = V * diag(sqrt(diag( D ))) * Y_est(:,:,k);
    end
    
end

if par.show_info
   fprintf('---------------omp ends,time cost:%.2f---------------\n',toc)
end
end