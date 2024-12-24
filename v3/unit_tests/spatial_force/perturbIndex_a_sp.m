function out = perturbIndex_a_sp(model,func,x,index_jj, index_ii)
% jt_index is the joint j w.r.t which the partial is taken
  step = 1e-20;
  m = length(x);
  ni = model.nv(index_ii); nj = model.nv(index_jj);
   
     for ind = 1:model.nv(index_jj)
           e = zeros(m,1);
           e(model.vinds{index_jj}(ind)) = 1;
           out_cs{ind} = func(x+1i*e*step);
     end  

     
     out.dfi_da_cs = zeros(6,nj);
     out.dfci_da_cs = zeros(6,nj);

    for ind=1:nj
        out.dfi_da_cs(:,ind) = imag(out_cs{ind}.f{index_ii})/step;
        out.dfci_da_cs(:,ind) = imag(out_cs{ind}.fC{index_ii})/step;
        
    end
end

