function X = vi_update(Yn,MIMO_info,X,Dic)

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


% update lambda: sparsity controller (precision) for the path loss
alpha_lam = K + X.alpha_lambda;
beta_lam = X.beta_lambda;
for k = 1:K
    beta_lam = beta_lam + real(X.x(:,k).*conj(X.x(:,k)) + diag(X.Var(:,:,k)));
end
X.lambda = alpha_lam./beta_lam;

% update xi: noise precision
alpha_xi = P*Q*K + X.alpha_xi;
beta_xi = X.beta_xi;
for k = 1:K
    ykn = reshape(Yn(:,:,k),P*Q,1);
    Dk_select = Dic.D(:,:,k);
    beta_xi = beta_xi + real(ykn'*ykn) + sum(   real(   diag(   ...
        Dk_select'*Dk_select*(   X.Var(:,:,k)...
        +X.x(:,k)*X.x(:,k)'   )   )   )   ) ...
        - 2*real(   ykn'*Dk_select*X.x(:,k)   );
end
X.xi = alpha_xi./beta_xi;

for k = 1:K
    Dk_select = Dic.D(:,:,k);
    lambda_select = X.lambda;
 
    % update X.x: the augmented path loss
    X.Var(:,:,k) = inv_for_XhX_plus_diag(   X.xi,Dk_select,lambda_select   ); % inv(X.xi*Dk'*Dk + diag(lambda))
    X.x(:,k) = X.xi * X.Var(:,:,k) * Dk_select' * reshape(Yn(:,:,k),P*Q,1);
end

end