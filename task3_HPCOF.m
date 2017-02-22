function [exp_data,s_HPCOF ,s_HPCOF_uptake,rho_HPCOF_a,pvalue_HPCOF_a,sse_HPCOF_a] =task3_HPCOF(model,fluxdata_file,fluxlist2retrieve_file,...
    flux2scale,iProblem,v_HPCOF,uptakeFlux,iUptake)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fluxdata = importdata(fluxdata_file);
exp_data = fluxdata.data(:,1);
%------------------------------------
% standard deviation of flux data
if size(fluxdata.data,2)==2
    exp_data_sd = fluxdata.data(:,2);
else
    exp_data_sd = [];
end
%------------------------------------

% initialize
fluxname    = fluxdata.textdata;
fluxname0   = importdata(fluxlist2retrieve_file);
s_HPCOF0     = zeros(size(fluxname0,1),1);
for i = 1:size(fluxname0,1)
    ii = find(strcmp(fluxname0{i},model.rxns));
    s_HPCOF0(i) = flux2scale*v_HPCOF(ii);
end
s_HPCOF = matchFlux(s_HPCOF0,iProblem);
function r = matchFlux(r0,iProblem)

switch iProblem
    
    case {11,12,13,14,15}
        r = r0;
        r(1) = abs(r(1)); % glucose uptake 
       

end
end
exp_data_uptake = uptakeFlux;
s_HOCOF_uptake = flux2scale*abs(v_HPCOF(iUptake));
% biomass yield
% v_HPCOF_BM = v_HPCOF(model.c==1);
% rate= v_HPCOF_BM*exp_data_uptake;
%% Correlation coefficient(rho), sum of squared error (SSE)

% Sum of seuqred error (SSE)
sse_ = @(f,y)(sum((y-f).^2))/norm(y);
sse_2=@(f,y)(norm(y-f))/norm(y);
% incorporate uptake flux into flux vector
exp_data_a = [exp_data_uptake;exp_data];
s_HPCOF_a=[s_HPCOF_uptake;s_HPCOF];
% x=fluxname0';
% y=[exp_data;s_HPCOF];
% bar(x,y);
% title('实验值和预测值对比');
% xlabel('反应名称');
% ylabel('代谢通量');
% performance measures of e-fmin 
[R P] = corrcoef(s_HOCOF_a,exp_data_a);
rho_HOCOF_a  = R(1,2); pvalue_HOCOF_a = P(1,2);
sse_HOCOF_a = sse_2(s_HOCOF_a,exp_data_a);
end

