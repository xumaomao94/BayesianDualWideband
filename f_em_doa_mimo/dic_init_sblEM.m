function Dic = dic_init_sblEM(MIMO_info,dense_AOD,dense_AOA,index_select)

%% System information
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

dense_grid = dense_AOA*dense_AOD;

dense_index_select = length(index_select);
index_select_wrt_AOD = ceil(index_select/dense_AOA); % [dense_index_select,1], pair with index_select_wrt_AOA
index_select_wrt_AOA = index_select - (index_select_wrt_AOD-1).*dense_AOA; % [dense_index_select,1]

%% Build the grid for aod/aoa
Dic.grid_AOD_eq = linspace(   -0.5 + 1/dense_AOD/2, 0.5 - 1/dense_AOD/2, dense_AOD   )'; % [dense_aod,1] points in [-0.5,0.5], since eqAOD = sin(AOD)/2
Dic.grid_AOA_eq = linspace(   -0.5 + 1/dense_AOA/2, 0.5 - 1/dense_AOA/2, dense_AOA   )';

%% Dictionary for vec(Yn)
Dic.D = zeros(Q*P,dense_index_select,K); % keep only [dense_index_select] among [dense_theta*dense_phi] columns
Dic.D_partial_AOD = zeros(Q*P,dense_index_select,K);
Dic.D_partial_AOA = zeros(Q*P,dense_index_select,K);

for k = 1:K
    f_k = K_select(k) * f_s / K_0;
    dic_Abs = exp(   -1i * 2 * pi * (1+f_k/f_c) * (0:Nt-1)' * Dic.grid_AOD_eq'   );
    deriv_dic_Abs = diag(   -1i * 2 * pi * (1+f_k/f_c) * (0:Nt-1)'   ) ... % derivative towards the equivalent AOD
                                        * dic_Abs;
	dic_Abs = F.'*dic_Abs;
    deriv_dic_Abs = F.'*deriv_dic_Abs;
                           
    dic_Ams = exp(   -1i * 2 * pi * (1+f_k/f_c) * (0:Nr-1)' * Dic.grid_AOA_eq'   );
    deriv_dic_Ams = diag(   -1i * 2 * pi * (1+f_k/f_c) * (0:Nr-1)'   ) ...
                                        * dic_Ams;
	dic_Ams = W.'*dic_Ams;
    deriv_dic_Ams = W.'*deriv_dic_Ams;
    
    if dense_grid == dense_index_select
        Dic.D(:,:,k) = kron(   dic_Abs, dic_Ams   );
        Dic.D_partial_AOD(:,:,k) = kron(   deriv_dic_Abs, dic_Ams   );
        Dic.D_partial_AOA(:,:,k) = kron(   dic_Abs, deriv_dic_Ams   );
    else
        Dic.D(:,:,k) = khatrirao(   dic_Abs(:,index_select_wrt_AOD), dic_Ams(:,index_select_wrt_AOA)   );
        Dic.D_partial_AOD(:,:,k) = khatrirao(   deriv_dic_Abs(:,index_select_wrt_AOD), dic_Ams(:,index_select_wrt_AOA)   );
        Dic.D_partial_AOA(:,:,k) = khatrirao(   dic_Abs(:,index_select_wrt_AOD), deriv_dic_Ams(:,index_select_wrt_AOA)   );
    end
end

Dic.bound_AOD = zeros(dense_AOD,2);
Dic.bound_AOA = zeros(dense_AOA,2);
if dense_grid == dense_index_select
    Dic.bound_AOD = ones(dense_AOD,1) * [-1/2/dense_AOD, 1/2/dense_AOD];
    Dic.bound_AOA = ones(dense_AOA,1) * [-1/2/dense_AOA, 1/2/dense_AOA];
else
    AOD_index_select = unique(index_select_wrt_AOD);
    AODeq_select = Dic.grid_AOD_eq(AOD_index_select);
    Dic.bound_AOD(AOD_index_select,:) = [   [-0.5-AODeq_select(1);... % [dense_AOD_select,2]
    (AODeq_select(1:end-1)-AODeq_select(2:end))/2],...
    [(AODeq_select(2:end)-AODeq_select(1:end-1))/2; 0.5-AODeq_select(end)]   ];
%     Dic.bound_AOD = max(Dic.bound_AOD,-1/2/dense_AOD); % in case of unstable update in the frist iteration
%     Dic.bound_AOD = min(Dic.bound_AOD,1/2/dense_AOD);

    AOA_index_select = unique(index_select_wrt_AOA);
    AOAeq_select = Dic.grid_AOA_eq(AOA_index_select);
    Dic.bound_AOA(AOA_index_select,:) = [   [-0.5-AOAeq_select(1);...
    (AOAeq_select(1:end-1)-AOAeq_select(2:end))/2],...
    [(AOAeq_select(2:end)-AOAeq_select(1:end-1))/2; 0.5-AOAeq_select(end)]   ];
%     Dic.bound_AOA = max(Dic.bound_AOA,-1/2/dense_AOA);
%     Dic.bound_AOA = min(Dic.bound_AOA,1/2/dense_AOA);
end


end