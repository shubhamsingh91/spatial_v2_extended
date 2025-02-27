function [glob,bod] = perturbIndex_v_sp(model,func,x,index_jj, index_ii)
% jt_index is the joint j w.r.t which the partial is taken
  step = 1e-20;
  m = length(x);
  ni = model.nv(index_ii); nj = model.nv(index_jj);
   
     for ind = 1:model.nv(index_jj)
           e = zeros(m,1);
           e(model.vinds{index_jj}(ind)) = 1;
           [glob_cs{ind},bod_cs{ind}] = func(x+1i*e*step);
     end  
     glob.dSdi_dvj_cs = zeros(6,ni,nj);
     glob.dpsidi_dvj_cs = zeros(6,ni,nj);
     glob.dBCi_dvj_cs = zeros(6,6,nj);
     glob.dfi_dqv_cs = zeros(6,nj);
     glob.dfci_dqv_cs = zeros(6,nj);

     bod.dfci_dqv_cs = zeros(6,nj);

    for ind=1:nj
        % global coordinates
        glob.dSdi_dvj_cs(:,:,ind) = imag(glob_cs{ind}.Sd{index_ii})/step;
        glob.dpsidi_dvj_cs(:,:,ind) = imag(glob_cs{ind}.psid{index_ii})/step;
        glob.dBCi_dvj_cs(:,:,ind) = imag(glob_cs{ind}.BC{index_ii})/step;
        glob.dfi_dqv_cs(:,ind) = imag(glob_cs{ind}.f{index_ii})/step;
        glob.dfci_dqv_cs(:,ind) = imag(glob_cs{ind}.fC{index_ii})/step;
        
        % body coordinates
        bod.dfCi_dvj_cs(:,ind) = imag(bod_cs{ind}.fC{index_ii})/step;

    end
    
   
    
end

