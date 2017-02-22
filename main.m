function main
% Select the problem number (iProblem) 

iProblem = 14;

% iProblem 
% --------------------------------------------------------------------------

% problem setting
[model,RXDHPCOF,fluxdata_file, fluxlist2retrieve_file,...
 uptakeFlux,iUptake,rxnName2scaleRXD,rxnName2scaleFluxes,flux2scale] = task1_setProblem(iProblem);
% simulation
%v_HPCOF= task2_long( model,RXDHPCOF,rxnName2scaleFluxes);
v_HPCOF=task2a_long(model,RXDHPCOF,rxnName2scaleFluxes)
% post-processing
[exp_data,s_HPCOF ,s_HPCOF_uptake,rho_HPCOF_a,pvalue_HPCOF_a,sse_HPCOF_a] =task3_long(model,fluxdata_file,fluxlist2retrieve_file,...
    flux2scale,iProblem,v_HPCOF,uptakeFlux,iUptake)

