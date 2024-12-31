clc; clear all;

%% testing numerical values
qnum = [1,1,1].'; mnum = [1,2,3].';
vnum = [1,1,1].'; fnum = [1,2,3,4,5,6].';

s_ring_num  = s_ring(qnum,vnum)

s_ring_v_deriv_num = s_ring_v_deriv(qnum,vnum)

dsm_dq_num = dSm_dq(qnum,mnum)
% 
% dSTf_dq_num = dSTf_dq(qnum,fnum)
% 
