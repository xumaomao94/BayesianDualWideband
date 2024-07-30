function [theta_pair,phi_pair,index_pair] = anglepair(varargin)
% input should be columnwise pairs of theta/phi
% 1. pair [theta_list,phi_list] by minimizing the error between [theta_init,phi_init]
% 2. when no [theta_init,phi_init], can be seen as one 
%    round KNN with initial value of [theta_list(:,1),phi_list(:,1)]


if nargin == 3
    theta_list = varargin{1};
    phi_list = varargin{2};
    MIMO_info = varargin{3};
    
    K_0 = MIMO_info.K_0;
    K_select = MIMO_info.K_select;
    f_s = MIMO_info.f_s;
    f_c = MIMO_info.f_c;
    
    theta_mean = theta_list(:,1);
    phi_mean = phi_list(:,1);
    L = length(theta_mean);
    
    index_pair = zeros(size(theta_list));
    
    for k = 2:size(theta_list,2)
        theta_cmp = theta_list(:,k);
        theta_cmp_shift = theta_cmp + 1./(1+K_select(k)*f_s/K_0/f_c).*[-1,0,1];
        phi_cmp = phi_list(:,k);
        phi_cmp_shift = phi_cmp + 1./(1+K_select(k)*f_s/K_0/f_c).*[-1,0,1];
        
        for i = 1:L
            % take into condition like -0.49 and 0.49, which are actually
            % close to each other
            
            [theta_dist,index_theta_shift] = min(   abs(theta_mean(i) - theta_cmp_shift),[],2   );
            [phi_dist,index_phi_shift] = min(   abs(phi_mean(i) - phi_cmp_shift),[],2   );
            
            dist = theta_dist.^2 + phi_dist.^2;
            
            [~,i_max] = min(dist); i_max = i_max(1);
            index_pair(i,k) = i_max;
            theta_mean(i) = theta_mean(i)*(k-1)/k + ...
                theta_cmp_shift(   i_max,index_theta_shift(i_max)   )/k;
            if theta_mean(i)<-0.5
                theta_mean(i) = theta_mean(i)+1;
            elseif theta_mean(i)>0.5
                theta_mean(i) = theta_mean(i)-1;
            end
            
            phi_mean(i) = phi_mean(i)*(k-1)/k + ...
                phi_cmp_shift(   i_max,index_phi_shift(i_max)   )/k;
            if phi_mean(i)<-0.5
                phi_mean(i) = phi_mean(i)+1;
            elseif phi_mean(i)>0.5
                phi_mean(i) = phi_mean(i)-1;
            end
        
            theta_cmp_shift(i_max,:) = 100; % take a large value so it won't be chosen next time
            phi_cmp_shift(i_max,:) = 100;
        end
        
    end
    
    theta_pair = theta_mean';
    phi_pair = phi_mean';
    
elseif nargin == 6
    theta_std = varargin{1};
    phi_std = varargin{2};
    theta_list = varargin{3};
    phi_list = varargin{4};
    alpha_list = varargin{5};
    MIMO_info = varargin{6};
    
    K = MIMO_info.K;
    K_0 = MIMO_info.K_0;
    K_select = MIMO_info.K_select;
    f_s = MIMO_info.f_s;
    f_c = MIMO_info.f_c;
    
    if size(theta_list,2) == K
        multiple_estimation = true;
    elseif size(theta_list,2) == 1
        multiple_estimation = false;
        K = 1;
    end
    
    L = length(theta_std);
    theta_pair = zeros(L,K);
    phi_pair = zeros(L,K);
    index_pair = zeros(L,K);
    
    for k = 1:K
        [~,index_select] = sort(   -abs(alpha_list(:,k))   );
        index_select = index_select(   1:min(L,length(index_select))   );
        
    	theta_cmp = theta_list(:,k);
        if L > length(theta_cmp)
            theta_cmp = [theta_cmp;zeros(L-length(theta_cmp),1)];
        end
        phi_cmp = phi_list(:,k);
        if L > length(phi_cmp)
            phi_cmp = [phi_cmp;zeros(L-length(phi_cmp),1)];
        end
        
        if multiple_estimation
            theta_cmp_shift = theta_cmp + [-1,0,1];
            phi_cmp_shift = phi_cmp + [-1,0,1];
        else
            theta_cmp_shift = theta_cmp + 1./(1+K_select(k)*f_s/K_0/f_c).*[-1,0,1];
            phi_cmp_shift = phi_cmp + 1./(1+K_select(k)*f_s/K_0/f_c).*[-1,0,1];
        end
        for i = 1:L
            [theta_dist,index_theta_shift] = min(   abs(theta_std(i) - theta_cmp_shift),[],2   );
            [phi_dist,index_phi_shift] = min(   abs(phi_std(i) - phi_cmp_shift),[],2   );
            
            dist = theta_dist.^2 + phi_dist.^2;
            [~,i_max] = min(dist); i_max = i_max(1);
            index_pair(i,k) = i_max;
            theta_pair(i,k) = theta_cmp_shift(   i_max,index_theta_shift(i_max)   );
            phi_pair(i,k) = phi_cmp_shift(   i_max,index_phi_shift(i_max)   );
            
            theta_cmp_shift(i_max,:) = 100; % take a large value so it won't be chosen next time
            phi_cmp_shift(i_max,:) = 100;
        end
    end
    

end