function [model,RXDHPCOF,fluxdata_file, fluxlist2retrieve_file,...
    uptakeFlux,iUptake,rxnName2scaleRXD,rxnName2scaleFluxes,flux2scale] = task1_setProblem(iProblem)

% read data

switch iProblem
    
    case {11,12}
        
        % network (mat file)
        modelName = 'yeast_5.21_MCISB';
        % uptake flux name
        uptakeFluxName = 'r_1714'; % D-glucose exchange
        % reactoin name to scale RXD
        rxnName2scaleRXD = 'r_1166'; % glucose transport 
        % settings for postprocessing
        fluxlist2retrieve_file = 'Lee_fluxname_model.txt';
        rxnName2scaleFluxes = uptakeFluxName;
        
        % datafiles
        if iProblem == 11
            genedata_file = 'Lee_75_genedata.txt';
            fluxdata_file = 'Lee_75_fluxdata.txt';
            uptakeFlux = 16.5;
        elseif iProblem == 12
            genedata_file = 'Lee_85_genedata.txt';
            fluxdata_file = 'Lee_85_fluxdata.txt';
            uptakeFlux = 11;
        end
        
        flux2scale = uptakeFlux;        
        load([modelName,'.mat'])
        
    case {13,14,15
            }

        % network (mat file)
        modelName = 'Ecoli_iAF1260';
        % uptake flux name
        uptakeFluxName = 'EX_glc(e)'; % D-glucose exchange
        % reactoin name to scale RXD
        rxnName2scaleRXD   = 'GLCptspp'; % D-glucose transport via PEP:Pyr PTS (periplasm)
        fluxlist2retrieve_file = 'reference_fluxname.txt';
        rxnName2scaleFluxes = uptakeFluxName; 
        if iProblem==13
        % settings for postprocessing
            genedata_file='data_reference.xlsx';
            fluxdata_file = 'reference_flux.txt';
            uptakeFlux = 9.2; 
        elseif iProblem==14
            genedata_file='data_NOX.xlsx';
            fluxdata_file = 'NOX_flux.txt';
            uptakeFlux = 11.7; 
        elseif iProblem==15
            genedata_file='data_ATP.xlsx';
            fluxdata_file = 'ATP_flux.txt';
            uptakeFlux = 15.6; 
        end
        flux2scale = uptakeFlux;        
        load([modelName,'.mat'])   
end

% find uptake flux
iUptake = find(strcmp(uptakeFluxName,model.rxns));
% % load transcriptomic data
    [genedata,genename] = xlsread(genedata_file);
    genename=genename(:,1);
     genename(1)=[];
      GXD_ave	= genedata(:,1);
%     genedata = importdata(genedata_file);
%    genename = genedata.textdata(:,1);
%     genename(1) = [];
%     GXD_ave	= genedata.data(:,1);
% x=GXD_ave;
% n=length(GXD_ave);
% r=10*rand([n 1]);
% y=x+0.1*(r-x);
 %GXD_ave=y;
      if size(genedata,2)==2
         GXD_sd = genedata(:,2);
     else
      GXD_sd = zeros(size(genedata,1),1);
     end
%     if size(genedata.data,2)==2
%         GXD_sd = genedata.data(:,2);
%    else
%        GXD_sd = zeros(size(genedata.data,1),1);
%     end
%% map GXD (gene expression data) to RXD (reaction expression data)
% The following functions and m-files were taken from Lee et al. and 
% modiifed by Hyun-Seob Song 
%   'convertGXD2RXD' 
%   'addGeneData'
%   'task1a_AandB.m'
%   'task1b_AorB_max.m'
%   'task1c_AorB_sum.m'

disp('Mapping gene expressions onto reactions ...')

% 'A or B' -> max(A,B) 
    [RXDmax_ave,RXDmax_sd] = convertGXD2RXD(model,genename,GXD_ave,GXD_sd);
    % sds 0 -> small
    if max(RXDmax_sd) ~= 0    
        RXDmax_sd(RXDmax_sd == 0) = min(RXDmax_sd(RXDmax_sd>0))/2;
    end
    % scale gene expression data (by rxnName2scaleRXD)
%       iRXDref = find(strcmp(rxnName2scaleRXD,model.rxns));
%     RXDmax_sd = RXDmax_sd/RXDmax_ave(iRXDref);
%     RXDmax_ave = RXDmax_ave/RXDmax_ave(iRXDref);

disp('... done')
disp(' ')
RXDHPCOF.ave = RXDmax_ave; RXDHPCOF.sd = RXDmax_sd;

%--------------------------------------------------------------------------
function [r,r_sd] = convertGXD2RXD(m,g,t,t_sd)
% kieran: 16 sep 11

r       = zeros(size(m.rxns));
r_sd    = zeros(size(m.rxns));

for k = 1:length(g)
    g{k} = strrep(g{k},'-','_');
end

for k = 1:length(m.rxns)
    
    ga = m.grRules{k};
    ga = strrep(ga,'-','_');
	ga = strrep(ga, '  ', ' ');
    ga = strrep(ga, '( ', '(');
    ga = strrep(ga, ' )', ')');
    ga = strrep(ga,'.','_');

    w = regexp(ga,'\<\w*\>','match'); 
    w = setdiff(w,{'and','or','AND','OR'});
    
    for kk = 1:length(w)
        
		j = find(strcmp(w{kk},g));
        if isempty(j) 
            n = 0; 
            if max(t_sd)==0
                n_sd = 0; 
            else
                n_sd = min(t_sd(t_sd>0))/2; 
            end
        else
            n = t(j);
            n_sd = t_sd(j);
            if length(n) >= 2 % if multiple values exist, take average 
                n = mean(n); 
                n_sd = norm(n_sd,2)/length(n_sd);
            end
        end
        %-------------------------
        if n <= 1.e-6, n = 0; 
        end; 
        ga = regexprep(ga,['\<',w{kk},'\>'],[num2str(n),'pm',num2str(n_sd)]);
    end
    
    [n,n_sd] = addGeneData(ga);
    r(k) = n;
    r_sd(k) = n_sd;
    
end
end

%--------------------------------------------------------------------------
function str = AandB(str1,str2) %#ok<DEFNU>

    %ApmB = '([0-9\.])+pm([0-9\.]+)';

    %FIX Daniel
    ApmB = '([0-9\.[\-\+]e])+pm([0-9\.[\-\+]e]+)';

    match_expr      = ApmB;
    m1              = eval(regexprep(str1,match_expr,'$1'));
    s1              = eval(regexprep(str1,match_expr,'$2'));
    m2              = eval(regexprep(str2,match_expr,'$1'));
    s2              = eval(regexprep(str2,match_expr,'$2'));

    [m,j] = min([m1,m2]);

    if j == 1
       s = s1;
    else
        s = s2;
   end
   
    str = [num2str(m),'pm',num2str(s)];

end

function str = AorB(str1,str2) %#ok<DEFNU>

    %ApmB = '([0-9\.])+pm([0-9\.]+)';

    %FIX Daniel
     ApmB = '([0-9\.[\-\+]e])+pm([0-9\.[\-\+]e]+)';

   match_expr = ApmB;
   m1= eval(regexprep(str1,match_expr,'$1'));
   s1= eval(regexprep(str1,match_expr,'$2'));
   m2= eval(regexprep(str2,match_expr,'$1'));
   s2= eval(regexprep(str2,match_expr,'$2'));

 [m,j] = max([m1,m2]);
 
 if j == 1
     s = s1;
 else
     s = s2;
 end
%  m = m1 + m2;
%  s = sqrt(s1^2 + s2^2);
str = [num2str(m),'pm',num2str(s)];
end

function [n,n_sd] = addGeneData(g)

% kieran: 22 july 11

n = nan;
n_sd = nan;
ApmB = '[0-9\.]+pm[0-9\.]+';
tries =0;
f_and = @(a,b)AandB(a,b);
f_or = @(a,b)AorB(a,b);
sizeG = size(g,2); 
if sizeG >= 25000
    disp('grRule is too long...')
end

if ~isempty(g)
    
    while isnan(n)
        tries = tries+1;
            if tries > 1000
                fprintf(1, 'warning: stuck at loop evaluating %s\n', g);
                break
            end
        
        try 
            
            match_expr      = '([0-9\.])+pm([0-9\.]+)';
            g_av            = regexprep(g,match_expr,'$1');
            g_sd            = regexprep(g,match_expr,'$2');
            n               = eval(g_av);
            n_sd            = eval(g_sd);
            
        catch %#ok<CTCH>
            
            % replace brackets
            match_expr = ['\((',ApmB,')\)'];
            replace_expr = '$1';
            g = regexprep(g,match_expr,replace_expr);
           
            % replace 'and'
            match_expr      = ['(',ApmB,') and (',ApmB,')'];
            replace_expr    = '${f_and($1,$2)}';
            g = regexprep(g,match_expr,replace_expr,'once');
            
            % replace 'or'
            match_expr      = ['(',ApmB,') or (',ApmB,')'];
            replace_expr    = '${f_or($1,$2)}';
            g = regexprep(g,match_expr,replace_expr,'once'); 
        end
    end
end
end
end



