function X = delayest(Yn,MIMO_info,X,Dic,varargin)

%% read the parameters
p = inputParser;
addRequired(p,'Yn',@isnumeric);
addRequired(p,'MIMO_info',@isstruct);
addRequired(p,'X',@isstruct);
addRequired(p,'Dic',@isstruct);
addOptional(p,'method',"search",@isstring); %"ls": least square, only for equally spaced K_select
addOptional(p,'dense_tau',1e4,@isnumeric); % useful if method == "search"
parse(p,Yn,MIMO_info,X,Dic,varargin{:});
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

dense_tau = par.dense_tau;

%% delay estimation
index_select = X.index_select;
dense_index_select = length(index_select);

X.delay = zeros(dense_index_select,1);
X.alpha = zeros(dense_index_select,1);

delayT = K_0/f_s;
tau_grid = linspace(0,delayT,dense_tau);

for ell = 1:dense_index_select
    aug_alpha = X.x(ell,:).';
    if strcmp(par.method,"search")
        D_tau = exp(   -1i*2*pi * f_s/K_0 * K_select' * tau_grid   ); % K * dense_tau
        [~,imax] = max(   abs(   D_tau'*aug_alpha   )   );
        X.delay(ell) = tau_grid(imax);
    elseif strcmp(par.method,"ls")
        X.delay(ell) = angle(   aug_alpha(1:end-1)'*aug_alpha(2:end)   )...
                        /(   - 2*pi * f_s/K_0 * (K_select(2)-K_select(1))   );
        X.delay(ell) = mod(X.delay(ell),delayT);
    end
    
    phase_shift_wrt_k = exp(   -1i*2*pi * f_s/K_0 * X.delay(ell) * K_select'    );
%     X.alpha(ell) = (phase_shift_wrt_k'*X.x(ell,:).')/(phase_shift_wrt_k'*phase_shift_wrt_k);
    X.alpha(ell) = (phase_shift_wrt_k'*X.x(ell,:).')/length(phase_shift_wrt_k);
    X.x(ell,:) = X.alpha(ell) * phase_shift_wrt_k.';
end

end