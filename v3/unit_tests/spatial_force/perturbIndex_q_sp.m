function [glob,bod] = perturbIndex_q_sp(model,func,x,index_jj, index_ii)
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
  
      % body coordinates
     bod.dSi_dqj_cs = zeros(6,ni,nj);
     bod.dSdi_dqj_cs = zeros(6,ni,nj);
     bod.dai_dqj_cs = zeros(6,nj);
     bod.dvJxsi_dqj_cs = zeros(6,ni,nj);
     bod.dpsidi_dqj_cs = zeros(6,ni,nj);
       bod.dIi_dqj_cs = zeros(6,6,nj);
       bod.dICi_dqj_cs = zeros(6,6,nj);
         bod.dai_dqj_cs = zeros(6,nj);
        bod.dIai_dqj_cs = zeros(6,nj);
   bod.dpsiddi_dqj_cs = zeros(6,ni,nj);
      bod.dBCi_dqj_cs = zeros(6,6,nj);
       bod.dfi_dqj_cs = zeros(6,nj);
       bod.dfCi_dqj_cs = zeros(6,nj);
      bod.dSiT_dqj_cs = zeros(ni,6,nj);
      bod.dvi_dqj_cs = zeros(6,nj);
      
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
    
    glob.diX0_dqj_cs =  imag(glob_cs{ind}.iXO{index_ii})/step;
    
    % body coordinates
    bod.dSi_dqj_cs(:,:,ind) = imag(bod_cs{ind}.S{index_ii})/step;
    bod.dSdi_dqj_cs(:,:,ind) = imag(bod_cs{ind}.Sd{index_ii})/step;
    bod.dvi_dqj_cs = imag(bod_cs{ind}.v{index_ii})/step;
    bod.dvJxsi_dqj_cs(:,:,ind) = imag(bod_cs{ind}.vJxS{index_ii})/step;
    bod.dpsidi_dqj_cs(:,:,ind) = imag(bod_cs{ind}.psid{index_ii})/step;
    bod.dIi_dqj_cs(:,:,ind) = imag(bod_cs{ind}.I{index_ii})/step;
    bod.dICi_dqj_cs(:,:,ind) = imag(bod_cs{ind}.IC{index_ii})/step;
    bod.dai_dqj_cs(:,ind) = imag(bod_cs{ind}.a{index_ii})/step;
    bod.dIai_dqj_cs(:,ind) = imag(bod_cs{ind}.Ia{index_ii})/step;
    bod.dpsiddi_dqj_cs(:,:,ind) = imag(bod_cs{ind}.psidd{index_ii})/step;
%     bod.dBCi_dqj_cs(:,:,ind) = imag(bod_cs{ind}.BC{index_ii})/step;
    bod.dfi_dqj_cs(:,ind) = imag(bod_cs{ind}.f{index_ii})/step;
%     bod.dfCi_dqj_cs(:,ind) = imag(bod_cs{ind}.fC{index_ii})/step;
    bod.dSiT_dqj_cs(:,:,ind) = imag(bod_cs{ind}.S{index_ii}.')/step;
%     bod.dvi_dqj_cs = imag(bod_cs{ind}.v{index_ii})/step;

    
    
end
end

