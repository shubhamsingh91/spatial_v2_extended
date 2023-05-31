function [glob] = perturbIndex_q(model,func,x,index_jj, index_ii,p)
% jt_index is the joint j w.r.t which the partial is taken
  step = 1e-20;
  m = length(x);
  ni = model.nv(index_ii); nj = model.nv(index_jj);
   
 for ind = 1:model.nv(index_jj)
       e = zeros(m,1);
       e(model.vinds{index_jj}(ind)) = 1;
       [glob_cs{ind},bod_cs{ind}] = func(x+1i*e*step);
 end  
    
     % global coordinates
      glob.dSi_dqj_cs = zeros(6,ni,nj);
     glob.dSdi_dqj_cs = zeros(6,ni,nj);
   glob.dvJxsi_dqj_cs = zeros(6,ni,nj);
   glob.dpsidi_dqj_cs = zeros(6,ni,nj);
       glob.dIi_dqj_cs = zeros(6,6,nj);
      glob.dICi_dqj_cs = zeros(6,6,nj);
         glob.dai_dqj_cs = zeros(6,nj);
        glob.dIai_dqj_cs = zeros(6,nj);
   glob.dpsiddi_dqj_cs = zeros(6,ni,nj);
      glob.dBCi_dqj_cs = zeros(6,6,nj);
       glob.dfi_dqj_cs = zeros(6,nj);
       glob.dfCi_dqj_cs = zeros(6,nj);
      glob.dSiT_dqj_cs = zeros(ni,6,nj);
      glob.dvi_dqj_cs = zeros(6,nj);
    
      glob.dSi_dqjp_cs = zeros(6,ni);
      
      glob.dSi_dqjp_cs = zeros(6,ni);
      glob.dIvi_dqj_cs = zeros(6,nj);
      glob.dImi_dqj_cs = zeros(6,nj);
      glob.dvisf_dqj_cs = zeros(6,nj);
      glob.dSiTf_dqj_cs = zeros(ni,nj);
      glob.dxii_dqj_cs = zeros(6,nj);
      glob.dgammai_dqj_cs = zeros(6,nj);
      
   glob.dSi_dqjp_cs =imag(glob_cs{p}.S{index_ii})/step;  
   
for ind=1:nj
    % global coordinates
    glob.dSi_dqj_cs(:,:,ind) = imag(glob_cs{ind}.S{index_ii})/step;
    glob.dSdi_dqj_cs(:,:,ind) = imag(glob_cs{ind}.Sd{index_ii})/step;
    glob.dvJxsi_dqj_cs(:,:,ind) = imag(glob_cs{ind}.vJxS{index_ii})/step;
    glob.dpsidi_dqj_cs(:,:,ind) = imag(glob_cs{ind}.psid{index_ii})/step;
    glob.dIi_dqj_cs(:,:,ind) = imag(glob_cs{ind}.I{index_ii})/step;
    glob.dICi_dqj_cs(:,:,ind) = imag(glob_cs{ind}.IC{index_ii})/step;
    glob.dai_dqj_cs(:,ind) = imag(glob_cs{ind}.a{index_ii})/step;
    glob.dIai_dqj_cs(:,ind) = imag(glob_cs{ind}.Ia{index_ii})/step;
    glob.dpsiddi_dqj_cs(:,:,ind) = imag(glob_cs{ind}.psidd{index_ii})/step;
    glob.dBCi_dqj_cs(:,:,ind) = imag(glob_cs{ind}.BC{index_ii})/step;
    glob.dfi_dqj_cs(:,ind) = imag(glob_cs{ind}.f{index_ii})/step;
    glob.dfCi_dqj_cs(:,ind) = imag(glob_cs{ind}.fC{index_ii})/step;
    glob.dSiT_dqj_cs(:,:,ind) = imag(glob_cs{ind}.S{index_ii}.')/step;
    glob.dvi_dqj_cs = imag(glob_cs{ind}.v{index_ii})/step;
    
    glob.dIvi_dqj_cs(:,ind) = imag(glob_cs{ind}.Iv{index_ii})/step;
    glob.dImi_dqj_cs(:,ind) = imag(glob_cs{ind}.Im{index_ii})/step;
    glob.dvisf_dqj_cs(:,ind) = imag(glob_cs{ind}.vxsf{index_ii})/step;   
    glob.dSiTf_dqj_cs(:,ind) = imag(glob_cs{ind}.STf{index_ii})/step;   
    glob.dxii_dqj_cs(:,ind) = imag(glob_cs{ind}.xi{index_ii})/step;   
    glob.dgammai_dqj_cs(:,ind) = imag(glob_cs{ind}.gamma{index_ii})/step;   

end
end

