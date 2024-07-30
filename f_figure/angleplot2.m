function [h1,h2] = angleplot3(aoa,aod,alpha,varargin)
% ----------------------------------------------------
% input: 
%   aoa:
%       aoa(i) is the aoa of the i-th path
%   aod:
%       aod(i) is the aod of the i-th path
%   alpha
%       alpha(i) is the abs pathloss of the i-th path
% ----------------------------------------------------
% output:
% x-axis: aoa
% y-axis: aod
% color: abs(alpha)

%%  parse parameters
p = inputParser;
addRequired(p,'aoa',@isnumeric);
addRequired(p,'aod',@isnumeric);
addRequired(p,'alpha',@(x) isnumeric(x) && length(x) == length(aoa) && length(x) == length(aod));
addOptional(p,'AOA_start',-0.5,@isnumeric); % equivalent aoa: -sin(pi)/2
addOptional(p,'AOA_end',0.5,@isnumeric);
addOptional(p,'AOD_start',-0.5,@isnumeric);
addOptional(p,'AOD_end',0.5,@isnumeric);
addOptional(p,'th',0.05,@isnumeric);
addOptional(p,'color_small_value',"ks",@isstring);
addOptional(p,'color_large_value',"gs",@isstring);

parse(p,aoa,aod,alpha,varargin{:});
par = p.Results;

AOA_start = par.AOA_start;
AOA_end = par.AOA_end;
AOD_start = par.AOD_start;
AOD_end = par.AOD_end;
th = par.th;

%% draw the angles
if aoa(1) < AOA_start
    aoa(1) = aoa(1) + AOA_end - AOA_start;
end
if aoa(end) > AOA_end
    aoa(end) = aoa(end) - (AOA_end - AOA_start);
end
if aod(1) < AOD_start
    aod(1) = aod(1) + AOD_end - AOD_start;
end
if aod(end) > AOD_end 
    aod(end) = aod(end) - (AOD_end - AOD_start);
end

alpha = abs(alpha);
indices = find(alpha>th);
h1 = scatter(aoa,aod,[],alpha);
hold on;
h2 = scatter(aoa(indices),aod(indices),[],alpha(indices));
xlim([AOA_start,AOA_end])
ylim([AOD_start,AOD_end])
xlabel('AOA');
ylabel('AOD');


end