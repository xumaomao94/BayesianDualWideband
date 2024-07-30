function Dic_update = em_update(Yn,MIMO_info,X,Dic,y_current)


%% read the parameters
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

y_temp = zeros(P*Q,K);
for k = 1:K
    y_temp(:,k) = Dic.D(:,:,k) * X.x(:,k);
end
Dic_update = Dic;

grid_phi_eq = Dic.grid_AOD_eq; % [dense,1]
grid_theta_eq = Dic.grid_AOA_eq;
dense_phi = length(grid_phi_eq);
dense_theta = length(grid_theta_eq);
dense_grid = dense_phi*dense_theta;

index_select = X.index_select;
dense_index_select = length(index_select);

% find which phi/theta in the grid should be updated
index_select_to_array = zeros(dense_phi*dense_theta,1);
index_select_to_array(index_select) = 1;

phi_select = convsum(index_select_to_array,dense_theta,1); % each phi_i control (i-1)n_theta+1 : ntheta
index_phi_select = find(phi_select>0);
dense_phi_select = length(index_phi_select);
indicator_phi_select = zeros(dense_phi,1);
indicator_phi_select(   index_phi_select   ) = 1:dense_phi_select;

theta_select = blocksum(index_select_to_array,dense_theta,1); % each theta_i control i:n_theta:end
index_theta_select = find(theta_select>0);
dense_theta_select = length(index_theta_select);
indicator_theta_select = zeros(dense_theta,1);
indicator_theta_select(   index_theta_select   ) = 1:dense_theta_select;

%% off-grid EM update
K11 = 0; K22 = 0; K21 = 0; J1 = 0; J2 = 0;
for k = 1:K
    Dk = Dic.D(:,:,k); % [Q*P , GtGr]
    Dk_p_phi = Dic.D_partial_AOD(:,:,k); % [Q*P , GtGr], column-wise derivative wrt phi
    Dk_p_theta = Dic.D_partial_AOA(:,:,k); 
    ynk = reshape(Yn(:,:,k),[P*Q,1]);
    E_xx = X.x(:,k) * X.x(:,k)'+ X.Var(:,:,k);
    
    K11 = K11 + real(   (Dk_p_phi'*Dk_p_phi).*E_xx   ); % length(index_select) * length(index_select)
    K22 = K22 + real(   (Dk_p_theta'*Dk_p_theta).*E_xx   );
    K21 = K21 + real(   (Dk_p_theta'*Dk_p_phi).*E_xx   );
    
    J1 = J1 - real(   (ynk'*Dk_p_phi).*X.x(:,k).'   ) ... % 1 * length(index_select)
          + real(   sum(   Dk_p_phi.*(E_xx*Dk').',1   )   );
    J2 = J2 - real(   (ynk'*Dk_p_theta).*X.x(:,k).'   ) ...
          + real(   sum(   Dk_p_theta.*(E_xx*Dk').',1   )   );
end

Kn11_column_recover = zeros(dense_phi_select,dense_index_select);
Kn22_column_recover = zeros(dense_theta_select,dense_index_select);
Kn21_column_recover = zeros(dense_theta_select,dense_index_select);
iphi = ceil(index_select/dense_theta);
itheta = index_select - (iphi-1).*dense_theta;
for i = 1:dense_index_select
    i_phi_select = indicator_phi_select(   iphi(i)   );
    i_theta_select = indicator_theta_select(   itheta(i)   );
    Kn11_column_recover(   i_phi_select,:   ) = Kn11_column_recover(   i_phi_select,:   )...
        + K11(   i,:   );
    Kn22_column_recover(   i_theta_select,:   ) = Kn22_column_recover(   i_theta_select,:   )...
        + K22(   i,:   );
    Kn21_column_recover(   i_theta_select,:   ) = Kn21_column_recover(   i_theta_select,:   )...
        + K21(   i,:   );
end
Kn11 = zeros(dense_phi_select,dense_phi_select);
Kn22 = zeros(dense_theta_select,dense_theta_select);
Kn21 = zeros(dense_theta_select,dense_phi_select);
Jn1 = zeros(1,dense_phi_select);
Jn2 = zeros(1,dense_theta_select);
for i = 1:dense_index_select
    i_phi_select = indicator_phi_select(   iphi(i)   );
    i_theta_select = indicator_theta_select(   itheta(i)   );
    Kn11(   :,i_phi_select   ) = Kn11(   :,i_phi_select   ) + ...
        Kn11_column_recover(   :,i   );
    Kn22(   :,i_theta_select   ) = Kn22(   :,i_theta_select   ) + ...
        Kn22_column_recover(   :,i   );
    Kn21(   :,i_phi_select   ) = Kn21(   :,i_phi_select   ) + ...
        Kn21_column_recover(   :,i   );
    Jn1(   1,i_phi_select   ) = Jn1(   1,i_phi_select   ) +...
        J1(   1,i   );
    Jn2(   1,i_theta_select   ) = Jn2(   1,i_theta_select   ) + ...
        J2(   1,i   );

end
Kn = [Kn11,   Kn21';
      Kn21,   Kn22];
Jn = [Jn1,Jn2];

sol = - Kn \ Jn'; % sol: [delta_phi;delta_theta]

%% update dictionary parameters
delta_AOD_eq = zeros(dense_phi,1);
delta_AOD_eq(index_phi_select) = min(   Dic.bound_AOD(index_phi_select,2),...
    max(   Dic.bound_AOD(index_phi_select,1), sol( 1:length(index_phi_select) )   )   );
delta_AOA_eq = zeros(dense_theta,1);
delta_AOA_eq(index_theta_select) = min(   Dic.bound_AOA(index_theta_select,2),...
    max(   Dic.bound_AOA(index_theta_select,1), sol( length(index_phi_select)+1:end )   )   );

% grid update
Dic_update.grid_AOD_eq = Dic.grid_AOD_eq + delta_AOD_eq;
Dic_update.grid_AOA_eq = Dic.grid_AOA_eq + delta_AOA_eq;

%% update the dictionary
Dic_update.D = zeros(Q*P,dense_index_select,K);
Dic_update.D_partial_AOD = zeros(Q*P,dense_index_select,K);
Dic_update.D_partial_AOA = zeros(Q*P,dense_index_select,K);
for k = 1:K
    f_k = K_select(k) * f_s / K_0;
    
    dic_Abs = zeros(Nt,dense_phi);
    dic_Abs(:,index_phi_select) = exp(   -1i * 2 * pi * (1+f_k/f_c) * (0:Nt-1)' * Dic_update.grid_AOD_eq(index_phi_select)'   );
    deriv_dic_Abs = zeros(Nt,dense_phi);
    deriv_dic_Abs(:,index_phi_select) = diag(   -1i * 2 * pi * (1+f_k/f_c) * (0:Nt-1)'   ) ... % derivative towards the equivalent AOD
                                        * dic_Abs(:,index_phi_select);
                                                         
    dic_Ams = zeros(Nr,dense_theta);
    dic_Ams(:,index_theta_select) = exp(   -1i * 2 * pi * (1+f_k/f_c) * (0:Nr-1)' * Dic_update.grid_AOA_eq(index_theta_select)'   );
    deriv_dic_Ams = zeros(Nr,dense_theta);
    deriv_dic_Ams(:,index_theta_select) = diag(   -1i * 2 * pi * (1+f_k/f_c) * (0:Nr-1)'   ) ...
                                        * dic_Ams(:,index_theta_select);
    
    iphi = ceil(index_select/dense_theta);
    itheta = index_select - (iphi-1).*dense_theta;
    
    Dic_update.D(:,:,k) = khatrirao(   F.'*dic_Abs(:,iphi), W.'*dic_Ams(:,itheta)   );
    Dic_update.D_partial_AOD(:,:,k) = khatrirao(   F.'*deriv_dic_Abs(:,iphi), W.'*dic_Ams(:,itheta)   );
    Dic_update.D_partial_AOA(:,:,k) = khatrirao(   F.'*dic_Abs(:,iphi), W.'*deriv_dic_Ams(:,itheta)   );
end

%% backtracking
y_update = zeros(P*Q,K);
for k = 1:K
    y_update(:,k) = Dic_update.D(:,:,k) * X.x(:,k);
end
if sum((y_update(:)-Yn(:)).^2) > sum((y_current(:)-Yn(:)).^2) && Dic_update.trackNum < 3
    Dic_update.trackNum = Dic_update.trackNum + 1;
    Dic_update.bound_AOD = Dic_update.bound_AOD/2;
    Dic_update.bound_AOA = Dic_update.bound_AOA/2;
    Dic_update = em_update(Yn,MIMO_info,X,Dic_update,y_current);
end
if sum((y_update(:)-Yn(:)).^2) > sum((y_temp(:)-Yn(:)).^2) && Dic.trackNum ~= 0
    Dic_update = Dic;
end

end