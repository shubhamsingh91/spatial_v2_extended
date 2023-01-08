function out = perturbIndex_v(model,func,x,index_jj, index_ii)
% jt_index is the joint j w.r.t which the partial is taken
  step = 1e-20;
  m = length(x);
  ni = model.nv(index_ii); nj = model.nv(index_jj);
   
     for ind = 1:model.nv(index_jj)
           e = zeros(m,1);
           e(model.vinds{index_jj}(ind)) = 1;
           out_cs{ind} = func(x+1i*e*step);
     end  
     out.dSdi_dvj_cs = zeros(6,ni,nj);
     out.dpsidi_dvj_cs = zeros(6,ni,nj);
     out.dBCi_dvj_cs = zeros(6,6,nj);

    for ind=1:nj
        out.dSdi_dvj_cs(:,:,ind) = imag(out_cs{ind}.Sd{index_ii})/step;
        out.dpsidi_dvj_cs(:,:,ind) = imag(out_cs{ind}.psid{index_ii})/step;
        out.dBCi_dvj_cs(:,:,ind) = imag(out_cs{ind}.BC{index_ii})/step;

    end
end

