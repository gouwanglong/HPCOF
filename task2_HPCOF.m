function v_HOCOF= task2_HPCOF( model,RXD,rxnName2scaleFluxes)
%using our algorithm to calculate each flux
% initialize
ge=RXD.ave;
%ge=ge/max(ge)
ge=(ge-min(ge))/(max(ge)-min(ge));
A=model.S;
ge(isnan(ge))=1;
%lb=model.lb;
%ub=model.ub;
for i = 1:length(ge)
      if model.lb(i) ~= 0
         model.lb(i) = -ge(i);
      end
      if model.ub(i) ~= 0
          model.ub(i) = ge(i);
     end
end
%model.lb(model.c==1)=0.001;
lb=model.lb;
ub=model.ub;
b=ge;
n=length(model.rxns);
iFlux2scale = find(strcmp(rxnName2scaleFluxes,model.rxns));
cvx_begin
    variable  v(n)
    minimize (sum(huber(v-b))+norm(v,1))
    subject to
            A*v==0
            lb<=v<=ub
cvx_end
 v_HPCOF = v/abs(v(iFlux2scale))
end

