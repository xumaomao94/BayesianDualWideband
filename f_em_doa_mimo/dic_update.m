function Dic = dic_update(MIMO_info,Dic)

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

%% Dictionary update
for k = 1:K
    f_k = K_select(k) * f_s / K_0;
    dic_Abs = exp(   -1i * 2 * pi * (1+f_k/f_c) * (0:Nt-1)' * Dic.grid_AOD_eq'   );
    deriv_dic_Abs = diag(   -1i * 2 * pi * (1+f_k/f_c) * (0:Nt-1)'   ) ... % derivative towards the equivalent AOD
                                        * dic_Abs;
                                    
    dic_Ams = exp(   -1i * 2 * pi * (1+f_k/f_c) * (0:Nr-1)' * Dic.grid_AOA_eq'   );
    deriv_dic_Ams = diag(   -1i * 2 * pi * (1+f_k/f_c) * (0:Nr-1)'   ) ...
                                        * dic_Ams;
                                    
    Dic.D(:,:,k) = kron(   F.'*dic_Abs, W.'*dic_Ams   );
    Dic.D_partial_AOD(:,:,k) = kron(   F.'*deriv_dic_Abs, W.'*dic_Ams   );
    Dic.D_partial_AOA(:,:,k) = kron(   F.'*dic_Abs, W.'*deriv_dic_Ams   );
end


end