function [Channel_est,Y_est] = em_offgrid_dualwideband(Yn,MIMO_info,varargin)

%%  parse parameters
p = inputParser;
addRequired(p,'Yn',@isnumeric);
addRequired(p,'MIMO_info',@isstruct);
addOptional(p,'dense_phi',128,@isnumeric);
addOptional(p,'dense_theta',128,@isnumeric);
addOptional(p,'init_method',"omp",@isstring); %"omp": omp init; "random": random init (large init dimension)
addOptional(p,'max_iter',500,@isnumeric);
addOptional(p,'tol',1e-12,@isnumeric);
addOptional(p,'threshold_start',0,@isnumeric);
addOptional(p,'threshold_end',5e-2,@isnumeric);
addOptional(p,'gridshift_on',true,@islogical);
addOptional(p,'whitening_or_not',true,@islogical);
addOptional(p,'angle_combine',true,@islogical);
addOptional(p,'Th_angle_combine',1/64,@isnumeric);
addOptional(p,'show_info',false,@islogical);
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

dense_phi = par.dense_phi;
dense_theta = par.dense_theta;
dense_grid = dense_phi * dense_theta;

max_iter = par.max_iter;
tol = par.tol;
Th = zeros(max_iter,1);
Th(1:min(max_iter,10)) = linspace(   par.threshold_start,par.threshold_end,...
                                min(max_iter,10)   );
Th(11:end) = par.threshold_end;

if par.show_info
   fprintf('------SBL-EM begins------\n')
   disp(par)
   tic
end

%% Probabilistic modeling, hyperparameter init, so that the model is sparsity-promoting
if strcmp(par.init_method,"omp")
    [Channel_init,~] = omp_dualwideband(Yn,MIMO_info,...
                                            'dense_phi',dense_phi,...
                                            'dense_theta',dense_theta,...
                                            'grid_sample_method',"sinspace",...
                                            'max_path_number',10,...
                                            'whitening_or_not',par.whitening_or_not,...
                                            'show_info',false,...
                                            'show_info_during_process',false);
	path_exist_or_not =  find(   sum(Channel_init.path_exist_or_not,2) > 0   );
    if length(path_exist_or_not) > 128 && par.gridshift_on
        [~,X.index_select] = maxk(   sum(abs(Channel_init.alpha_aug),2),128   );
        X.index_select = sort(X.index_select);
    else
        X.index_select = path_exist_or_not; % keep a path as long as it exist in one subchannel
    end
    
elseif strcmp(par.init_method,"random")
    X.index_select = (1:dense_grid)';
end
dense_index_select = length(X.index_select);

% whitening
if par.whitening_or_not
    [V,D] = eig(MIMO_info.W.'*conj(MIMO_info.W));
    whiteningM = diag(1./sqrt(diag( D )))*V';
    Yn = whiteningM * reshape(Yn,MIMO_info.Q,[]);
    Yn = reshape(Yn,[MIMO_info.Q,MIMO_info.P,MIMO_info.K]);
    MIMO_info.W = MIMO_info.W * whiteningM.';
end

% noise precision, gamma
X.alpha_xi = 1e-6;
X.beta_xi = 1e-6;
X.xi = X.alpha_xi/X.beta_xi;
% precision for the path loss, gamma
X.alpha_lambda = 1e-6 * ones(dense_index_select,1);
X.beta_lambda = 1e-6 * ones(dense_index_select,1);
X.lambda = X.alpha_lambda./X.beta_lambda;
% path loss, gaussian with zero mean and precision lambda
if strcmp(par.init_method,"omp")
    X.x = Channel_init.alpha_aug(X.index_select,:);
elseif strcmp(par.init_method,"random")
    X.x = zeros(dense_index_select,1);
end

X.Var = zeros(dense_index_select,dense_index_select,K);
for k = 1:K
    X.Var(:,:,k) = eye(dense_index_select,dense_index_select);
end

%% Build the dictionary, and derivatives
Dic = dic_init_sblEM(MIMO_info,dense_phi,dense_theta,X.index_select);

%% Channel parameter estimation
y_est = zeros(P*Q,K);
rse = 0;
for i = 1:max_iter
    if par.show_info_during_process
        if mod(i,10) == 0 || i == 1
             fprintf('iteration:%d, rse:%.12f, number of path:%d, time cost:%.2f\n',...
                i, rse ,length(X.index_select), toc)
        end
    end
    
    X = vi_update(Yn,MIMO_info,X,Dic);

    if par.gridshift_on
        Dic.trackNum = 0;
        Dic = em_update(Yn,MIMO_info,X,Dic,y_est);
    end

    if par.angle_combine && par.gridshift_on
        [X,Dic] = gridprune(dense_phi,dense_theta,X,Dic,Th(i),max_iter,i,"AngleCombine_On",par.Th_angle_combine);
    else
        [X,Dic] = gridprune(dense_phi,dense_theta,X,Dic,Th(i),max_iter,i);
    end

    if par.show_info_during_process
        % draw only the first subcarrier
        if exist('h2')
            delete(h1)
            delete(h2)
        end
        iphi = ceil(X.index_select/dense_theta);
        itheta = X.index_select - (iphi-1).*dense_theta;
        [h1,h2] = angleplot3(   Dic.grid_AOA_eq(itheta), Dic.grid_AOD_eq(iphi),...
            mean(abs(X.x),2),...
            'color_small_value',"ko",...
            'color_large_value',"bo");
    end
    
    y_temp = y_est;
    for k = 1:K
        y_est(:,k) = Dic.D(:,:,k) * X.x(:,k);
    end
    dist = y_est(:) - y_temp(:);
    rse = dist'*dist/(y_est(:)'*y_est(:));
    if rse < tol
        if par.show_info_during_process
            fprintf('iteration:%d, rse:%.12f, number of path:%d, time cost:%.2f\n',...
                    i, rse ,length(X.index_select), toc)
        end
        break
    end
end

%% estimate the delay
X = delayest(Yn,MIMO_info,X,Dic,'dense_tau',1e6);

%% Rebuild the channel
iphi = ceil(X.index_select/dense_theta);
itheta = X.index_select - (iphi-1).*dense_theta;
if par.show_info
    if exist('h2')
        delete(h1)
        delete(h2)
    end
    fprintf('------SBL-EM converges, total time cost:%.2f------\n',toc)
    angleplot3(   Dic.grid_AOA_eq(itheta), Dic.grid_AOD_eq(iphi),...
        abs(X.alpha),...
        'color_small_value',"ko",...
        'color_large_value',"bo"   );
end

Channel_est.phi = Dic.grid_AOD_eq(iphi);
Channel_est.theta = Dic.grid_AOA_eq(itheta);
Channel_est.alpha = X.alpha;
Channel_est.tau = X.delay;
Channel_est.alpha_aug = X.x;

Channel_est.H = zeros(Nr,Nt,K);
Y_est = zeros(Q,P,K);
for k = 1:K
    f_k = K_select(k) * f_s/K_0;
    Channel_est.H(:,:,k) = exp(   -1i.*2.*pi .* (1+f_k/f_c) .*...
        (0:Nr-1)' * Channel_est.theta'   ) * diag(Channel_est.alpha_aug(:,k)) *...
        exp(   -1i.*2.*pi .* (1+f_k/f_c) .* (0:Nt-1)' * Channel_est.phi'   ).';
    Y_est(:,:,k) = MIMO_info.W.' * Channel_est.H(:,:,k) * F;
    if par.whitening_or_not
        Y_est(:,:,k) = V * diag(sqrt(diag( D ))) * Y_est(:,:,k);
    end
end

end