function [glob, bod]= perturbIndex_a_sp(model,func,x,index_jj, index_ii)
% jt_index is the joint j w.r.t which the partial is taken
  step = 1e-20;
  m = length(x);
  ni = model.nv(index_ii); nj = model.nv(index_jj);
   
     for ind = 1:model.nv(index_jj)
           e = zeros(m,1);
           e(model.vinds{index_jj}(ind)) = 1;
           [glob_cs{ind},bod_cs{ind}] = func(x+1i*e*step);
     end  

     
     glob.dfi_da_cs = zeros(6,nj);
     glob.dfci_da_cs = zeros(6,nj);

    for ind=1:nj
        % global coords
        glob.dfi_da_cs(:,ind) = imag(glob_cs{ind}.f{index_ii})/step;
        glob.dfci_da_cs(:,ind) = imag(glob_cs{ind}.fC{index_ii})/step;
    
        % body coords
        bod.dfci_da_cs(:,ind) = imag(bod_cs{ind}.fC{index_ii})/step;

    end
end

