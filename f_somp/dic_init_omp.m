function [Dic,Dic_sqr] = dic_init_omp(MIMO_info,dense_phi,dense_theta,k_carrier,grid_sample_method)

%% System Information
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

%% Build the grid (AOD and AOA)
if strcmp(grid_sample_method,"anglespace")
    grid_phi = linspace(   -pi/2+pi/dense_phi/2, pi/2-pi/dense_phi/2, dense_phi  );
    grid_phi_eq = sin(grid_phi)./2; % [1,dense_phi]
    grid_theta = linspace(   -pi/2+pi/dense_theta/2, pi/2-pi/dense_theta/2, dense_theta  );
    grid_theta_eq = sin(grid_theta)./2; % [1,dense_theta]
elseif strcmp(grid_sample_method,"sinspace")
    grid_phi_eq = linspace(   -0.5 + 1/dense_phi/2, 0.5 - 1/dense_phi/2, dense_phi   ); % dense_phi points in [-0.5,0.5], since eqAOD = sin(AOD)/2
    grid_theta_eq = linspace(   -0.5 + 1/dense_theta/2, 0.5 - 1/dense_theta/2, dense_theta   );
end

%% Build the dictionary
f_k = K_select(k_carrier) * f_s / K_0;
Dic_phi = F.' * exp(   -1i.*2.*pi .* (1+f_k/f_c) .* (0:Nt-1)' .* grid_phi_eq   ); % P * dense_phi
Dic_theta = W.' * exp(   -1i.*2.*pi .* (1+f_k/f_c) .* (0:Nr-1)' .* grid_theta_eq   ); % Q * dense_theta

Dic = kron(Dic_phi,Dic_theta);
Dic_sqr = kron(  sum( real(conj(Dic_phi).*Dic_phi),1),...
    sum( real(conj(Dic_theta).*Dic_theta),1)   )';
end