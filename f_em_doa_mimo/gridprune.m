function [X,Dic] = gridprune(dense_phi,dense_theta,X,Dic,th,max_iter,iter,varargin)

%% combine collapsed angles
index_select_wrt_phi = ceil(X.index_select/dense_theta); % [dense_index_select,1], pair with index_select_wrt_theta
index_select_wrt_theta = X.index_select - (index_select_wrt_phi-1).*dense_theta; % [dense_index_select,1]
index_select_wrt_phi_unique = unique(index_select_wrt_phi);
phi_dist = Dic.grid_AOD_eq(index_select_wrt_phi_unique(2:end)) - ...
    Dic.grid_AOD_eq(index_select_wrt_phi_unique(1:end-1));
index_phi_to_combine_unique = find(phi_dist<9.7656e-04);%1e-6);
for i = 1:length(index_phi_to_combine_unique)
    index_phi_to_combine = index_select_wrt_phi ==...
        index_select_wrt_phi_unique(index_phi_to_combine_unique(i));
    index_phi_to_copy = find(index_select_wrt_phi ==...
        index_select_wrt_phi_unique(index_phi_to_combine_unique(i)+1),1);
    index_to_copy = (index_select_wrt_phi(index_phi_to_copy)-1)...
        *dense_theta + index_select_wrt_theta(index_phi_to_combine);
    if ~ismember(index_to_copy,X.index_select)
        Dic.grid_AOD_eq(index_select_wrt_phi_unique(index_phi_to_combine_unique)) = 0;
        X.index_select(index_phi_to_combine) = index_to_copy;
    end
end
if isempty(index_phi_to_combine_unique)
    if Dic.grid_AOD_eq(index_select_wrt_phi_unique(1)) < -0.5
        Dic.grid_AOD_eq = [Dic.grid_AOD_eq(index_select_wrt_phi_unique(1)+1:end);...
            Dic.grid_AOD_eq(1:index_select_wrt_phi_unique(1)-1);...
            Dic.grid_AOD_eq(index_select_wrt_phi_unique(1))+1];
        X.index_select = X.index_select - dense_theta*index_select_wrt_phi_unique(1);
        index_phi_to_combine_unique = 1;
        index_phi_to_combine = index_select_wrt_phi ==...
            index_select_wrt_phi_unique(index_phi_to_combine_unique(1));
        X.index_select(index_phi_to_combine) = X.index_select(index_phi_to_combine)...
            + dense_theta*(dense_phi);
    elseif Dic.grid_AOD_eq(index_select_wrt_phi_unique(end)) > 0.5
        Dic.grid_AOD_eq = [Dic.grid_AOD_eq(index_select_wrt_phi_unique(end))-1;...
            Dic.grid_AOD_eq(index_select_wrt_phi_unique(end)+1:end);...
            Dic.grid_AOD_eq(1:index_select_wrt_phi_unique(end)-1)];
        X.index_select = X.index_select + dense_theta*(dense_phi-index_select_wrt_phi_unique(end)+1);
        index_phi_to_combine_unique = length(index_select_wrt_phi_unique);
        index_phi_to_combine = index_select_wrt_phi ==...
            index_select_wrt_phi_unique(index_phi_to_combine_unique(1));
        X.index_select(index_phi_to_combine) = X.index_select(index_phi_to_combine)...
            - dense_theta*dense_phi;
    end
end

index_select_wrt_phi = ceil(X.index_select/dense_theta); % [dense_index_select,1], pair with index_select_wrt_theta
index_select_wrt_theta = X.index_select - (index_select_wrt_phi-1).*dense_theta; % [dense_index_select,1]
index_select_wrt_theta_unique = unique(index_select_wrt_theta);
theta_dist = Dic.grid_AOA_eq(index_select_wrt_theta_unique(2:end)) - ...
    Dic.grid_AOA_eq(index_select_wrt_theta_unique(1:end-1));
index_theta_to_combine_unique = find(theta_dist<9.7656e-04);%1e-6);
for i = 1:length(index_theta_to_combine_unique)
    index_theta_to_combine = index_select_wrt_theta ==...
        index_select_wrt_theta_unique(index_theta_to_combine_unique(i));
    index_theta_to_copy = find(index_select_wrt_theta ==...
        index_select_wrt_theta_unique(index_theta_to_combine_unique(i)+1),1);
    index_to_copy = (index_select_wrt_phi(index_theta_to_combine)-1)...
        *dense_theta + index_select_wrt_theta(index_theta_to_copy);
    if ~ismember(index_to_copy,X.index_select)
        Dic.grid_AOA_eq(index_select_wrt_theta_unique(index_theta_to_combine_unique)) = 0;
        X.index_select(index_theta_to_combine) = index_to_copy;
    end
end
if isempty(index_theta_to_combine_unique)
    if Dic.grid_AOA_eq(index_select_wrt_theta_unique(1)) < -0.5
        Dic.grid_AOA_eq = [Dic.grid_AOA_eq(index_select_wrt_theta_unique(1)+1:end);...
            Dic.grid_AOA_eq(1:index_select_wrt_theta_unique(1)-1);...
            Dic.grid_AOA_eq(index_select_wrt_theta_unique(1))+1];
        X.index_select = X.index_select - index_select_wrt_theta_unique(1);
        index_theta_to_combine_unique = 1;
        index_theta_to_combine = index_select_wrt_theta ==...
            index_select_wrt_theta_unique(index_theta_to_combine_unique(1));
        X.index_select(index_theta_to_combine) = X.index_select(index_theta_to_combine)...
            + (dense_theta);
    elseif Dic.grid_AOA_eq(index_select_wrt_theta_unique(end)) > 0.5
        Dic.grid_AOA_eq = [Dic.grid_AOA_eq(index_select_wrt_theta_unique(end))-1;...
            Dic.grid_AOA_eq(index_select_wrt_theta_unique(end)+1:end);...
            Dic.grid_AOA_eq(1:index_select_wrt_theta_unique(end)-1)];
        X.index_select = X.index_select + (dense_theta-index_select_wrt_theta_unique(end)+1);
        index_theta_to_combine_unique = length(index_select_wrt_theta_unique);
        index_theta_to_combine = index_select_wrt_theta ==...
            index_select_wrt_theta_unique(index_theta_to_combine_unique(1));
        X.index_select(index_theta_to_combine) = X.index_select(index_theta_to_combine)...
            - dense_theta;
    end
end

if ~isempty(index_phi_to_combine_unique) || ~isempty(index_theta_to_combine_unique)
    [X.index_select,perm] = sort(X.index_select);
    X.x = X.x(perm,:);
    X.lambda = X.lambda(perm);
    X.Var = X.Var(perm,perm,:);
    Dic.D = Dic.D(:,perm,:);
    Dic.D_partial_AOD = Dic.D_partial_AOD(:,perm,:);
    Dic.D_partial_AOA = Dic.D_partial_AOA(:,perm,:);
end

%% pruning unimportant components
index_discard_existing_basis = find(   mean(   abs(X.x) ,2   ) <= th   );
index_discard_in_grid = X.index_select(   index_discard_existing_basis   );

%% whether combine nearby angles into one
if ~isempty(varargin) && strcmp(varargin{1},"AngleCombine_On")
    
    index_select_wrt_phi = ceil(X.index_select/dense_theta); % [dense_index_select,1], pair with index_select_wrt_theta
    index_select_wrt_theta = X.index_select - (index_select_wrt_phi-1).*dense_theta; % [dense_index_select,1]
    phi_dist = Dic.grid_AOD_eq(index_select_wrt_phi); % [dense_index_select,1]
    phi_dist = triu(   phi_dist - phi_dist'   );
    theta_dist = Dic.grid_AOA_eq(index_select_wrt_theta);
    theta_dist = triu(   theta_dist - theta_dist'   );
    
    distance = sqrt(   phi_dist.^2 + theta_dist.^2   );
    
    if ~isempty(varargin{2})
        Th_anglecombine = min(varargin{2},1/256) + min(iter,20)/32/min(max_iter,20);
        [index_to_combine1,index_to_combine2] = find(   distance>0 & distance<Th_anglecombine   );
    else
        Th_anglecombine = 1/256 + min(iter,20)/64/min(max_iter,20);
        [index_to_combine1,index_to_combine2] = find(   distance>0 & distance<Th_anglecombine   );
    end
    
    for i = 1:length(index_to_combine1)
        if mean(   abs(   X.x(index_to_combine1(i),:)   )   )...
                <= mean(   abs(   X.x(index_to_combine2(i),:)   )   )
            if ~ismember(   index_to_combine1(i),index_discard_existing_basis   ) && ...
                    ~ismember(   index_to_combine2(i),index_discard_existing_basis   )
                index_discard_existing_basis = [index_discard_existing_basis;...
                    index_to_combine1(i)];
                index_discard_in_grid = [index_discard_in_grid;...
                    X.index_select(   index_to_combine1(i)   )];
            end
        else
            if ~ismember(   index_to_combine1(i),index_discard_existing_basis   ) && ...
                    ~ismember(   index_to_combine2(i),index_discard_existing_basis   )
                index_discard_existing_basis = [index_discard_existing_basis;...
                    index_to_combine2(i)];
                index_discard_in_grid = [index_discard_in_grid;...
                    X.index_select(   index_to_combine2(i)   )];
            end
        end
    end
    
    %% select only one path from the same AOA/AOD
    if iter > 10
        if mod(iter,2) == 0
            index_unique_phi = unique(index_select_wrt_phi);
            for i = 1:length(index_unique_phi)
                index_same_phi = find(index_select_wrt_phi == index_unique_phi(i));
                if length(index_same_phi) > 1
                    index_theta_to_compare = index_select_wrt_theta(index_same_phi);
                    index_to_compare = (index_unique_phi(i)-1)*dense_theta + index_theta_to_compare;
                    index_to_compare_existing_basis...
                        = find(ismember(X.index_select,index_to_compare));
                    [~,index_to_discard_same_phi] = min(abs(X.x(index_to_compare_existing_basis)));
                    index_to_discard_same_phi...
                        = index_to_compare_existing_basis(index_to_discard_same_phi(1));
                    if ~ismember(   index_to_discard_same_phi,...
                            index_discard_existing_basis   )
                        index_discard_existing_basis = [index_discard_existing_basis;...
                            index_to_discard_same_phi(1)];
                        index_discard_in_grid = [index_discard_in_grid;...
                        X.index_select(   index_to_discard_same_phi   )];
                    end
                end
            end
        else
        index_unique_theta = unique(index_select_wrt_theta);
            for i = 1:length(index_unique_theta)
                index_same_theta = find(index_select_wrt_theta == index_unique_theta(i));
                if length(index_same_theta) > 1
                    index_phi_to_compare = index_select_wrt_phi(index_same_theta);
                    index_to_compare = (index_phi_to_compare-1)*dense_theta + index_unique_theta(i);
                    index_to_compare_existing_basis...
                        = find(ismember(X.index_select,index_to_compare));
                    [~,index_to_discard_same_theta] = min(abs(X.x(index_to_compare_existing_basis)));
                    index_to_discard_same_theta...
                        = index_to_compare_existing_basis(index_to_discard_same_theta(1));
                    if ~ismember(   index_to_discard_same_theta,...
                            index_discard_existing_basis   )
                        index_discard_existing_basis = [index_discard_existing_basis;...
                            index_to_discard_same_theta(1)];
                        index_discard_in_grid = [index_discard_in_grid;...
                        X.index_select(   index_to_discard_same_theta(1)   )];
                    end
                end
            end
        end
    end
    
end

X.index_select = setdiff(   X.index_select,index_discard_in_grid   );
dense_grid = dense_phi*dense_theta;
dense_index_select = length(X.index_select);

index_select_wrt_phi = ceil(X.index_select/dense_theta); % [dense_index_select,1], pair with index_select_wrt_theta
index_select_wrt_theta = X.index_select - (index_select_wrt_phi-1).*dense_theta; % [dense_index_select,1]
index_phi_select = unique(index_select_wrt_phi);
index_theta_select = unique(index_select_wrt_theta);

%% discard dictionary
Dic.D(:,index_discard_existing_basis,:) = [];
Dic.D_partial_AOD(:,index_discard_existing_basis,:) = [];
Dic.D_partial_AOA(:,index_discard_existing_basis,:) = [];
temp = zeros(dense_phi,1);
temp(index_phi_select) = Dic.grid_AOD_eq(index_phi_select);
Dic.grid_AOD_eq = temp;
temp = zeros(dense_theta,1);
temp(index_theta_select) = Dic.grid_AOA_eq(index_theta_select);
Dic.grid_AOA_eq = temp;

AODeq_select = Dic.grid_AOD_eq(index_phi_select);
Dic.bound_AOD = zeros(dense_phi,2);
Dic.bound_AOD(index_phi_select,:) = [   [(AODeq_select(end)-1-AODeq_select(1))/2;...
    (AODeq_select(1:end-1)-AODeq_select(2:end))/2],...
    [(AODeq_select(2:end)-AODeq_select(1:end-1))/2; (AODeq_select(1)+1-AODeq_select(end))/2]   ];
if iter >= 1
    Dic.bound_AOD = Dic.bound_AOD/(max_iter)*(max_iter+1-iter);
    Dic.bound_AOD = max(Dic.bound_AOD,-1/max(dense_phi,dense_theta)/2);
    Dic.bound_AOD = min(Dic.bound_AOD,1/max(dense_phi,dense_theta)/2);
end
    
AOAeq_select = Dic.grid_AOA_eq(index_theta_select);
Dic.bound_AOA = zeros(dense_theta,2);
Dic.bound_AOA(index_theta_select,:) = [   [(AOAeq_select(end)-1-AOAeq_select(1))/2;...
    (AOAeq_select(1:end-1)-AOAeq_select(2:end))/2],...
    [(AOAeq_select(2:end)-AOAeq_select(1:end-1))/2; (AOAeq_select(1)+1-AOAeq_select(end))/2]   ];
if iter >= 1
    Dic.bound_AOA = Dic.bound_AOA/(max_iter)*(max_iter+1-iter);
    Dic.bound_AOA = max(Dic.bound_AOA,-1/max(dense_phi,dense_theta)/2);
    Dic.bound_AOA = min(Dic.bound_AOA,1/max(dense_phi,dense_theta)/2);
end

%% discard model parameters
X.x(index_discard_existing_basis,:) = [];
X.Var(index_discard_existing_basis,:,:) = [];
X.Var(:,index_discard_existing_basis,:) = [];

X.alpha_lambda(index_discard_existing_basis) = [];
X.beta_lambda(index_discard_existing_basis) = [];
X.lambda(index_discard_existing_basis) = [];

end