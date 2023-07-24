function [robot] = Jac_build(robot, q, qd, qdd, lam, foot)

if foot == 1 % front-foot
    robot.J = jac_fr_ft(q(1), q(2), q(3), q(4), q(5), q(6), q(7));
    robot.Jd = jacd_fr_ft(q(1), q(2), q(3), q(4), q(5), q(6), q(7), qd(1), qd(2), qd(3), qd(4), qd(5), qd(6), qd(7));
    robot.JTlam_FO = JTlam_fr_ft_FO_q(q(1),q(2),q(3),q(4),q(5),q(6),q(7),lam(1),lam(2));
    robot.Jqdd_FO_q = Jqdd_fr_ft_FO_q(q(1),q(2),q(3),q(4),...
                    q(5),q(6),q(7),qdd(1),qdd(2),...
                    qdd(3),qdd(4),qdd(5),qdd(6),qdd(7));
    robot.Jdqd_FO_q = Jdqd_fr_ft_FO_q(q(1),q(2),q(3),q(4),q(5),q(6),q(7),...
                qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7));  
   robot.Jdqd_FO_v = Jdqd_fr_ft_FO_v(q(1),q(2),q(3),q(4),q(5),q(6),q(7),...
                 qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7));   
  
   robot.JTlam_SOq = JTlam_fr_ft_SO_q(q(1),q(2),q(3),q(4),q(5),q(6),q(7),lam(1),lam(2));
  
  robot.Jqdd_SO_q = Jqdd_fr_ft_SO_q(q(1),q(2),q(3),q(4),q(5),q(6),q(7),qdd(1),qdd(2),...
                    qdd(3),qdd(4),qdd(5),qdd(6),qdd(7));
                
  robot.Jdqd_SO_q = Jdqd_fr_ft_SO_q(q(1),q(2),q(3),q(4),q(5),q(6),q(7),...
                  qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7)); 
             
  robot.Jdqd_SO_v = Jdqd_fr_ft_SO_v(q(1),q(2),q(3),q(4),q(5),q(6),q(7),...
                  qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7));  

  robot.Jdqd_SO_vq = Jdqd_fr_ft_SO_vq(q(1),q(2),q(3),q(4),q(5),q(6),q(7),...
                  qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7));  
          
  robot.J_FO_q = J_front_ft_FO_q(q(1),q(2),q(3),q(4),q(5),q(6),q(7));
  robot.JT_FO_q = JT_front_ft_FO_q(q(1),q(2),q(3),q(4),q(5),q(6),q(7));

elseif foot == 2 % back-foot

    robot.J = jac_ba_ft(q(1), q(2), q(3), q(4), q(5), q(6), q(7));
    robot.Jd = jacd_b_ft(q(1), q(2), q(3), q(4), q(5), q(6), q(7), qd(1), qd(2), qd(3), qd(4), qd(5), qd(6), qd(7));
    robot.JTlam_FO = JTlam_b_ft_FO_q(q(1),q(2),q(3),q(4),q(5),q(6),q(7),lam(1),lam(2));
    robot.Jqdd_FO_q = Jqdd_b_ft_FO_q(q(1),q(2),q(3),q(4),...
                    q(5),q(6),q(7),qdd(1),qdd(2),...
                    qdd(3),qdd(4),qdd(5),qdd(6),qdd(7));
    robot.Jdqd_FO_q = Jdqd_b_ft_FO_q(q(1),q(2),q(3),q(4),q(5),q(6),q(7),...
                qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7));  
   robot.Jdqd_FO_v = Jdqd_b_ft_FO_v(q(1),q(2),q(3),q(4),q(5),q(6),q(7),...
                 qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7));   
             
   robot.JTlam_SOq = JTlam_b_ft_SO_q(q(1),q(2),q(3),q(4),q(5),q(6),q(7),lam(1),lam(2));
  
   robot.Jqdd_SO_q = Jqdd_b_ft_SO_q(q(1),q(2),q(3),q(4),q(5),q(6),q(7),qdd(1),qdd(2),...
                    qdd(3),qdd(4),qdd(5),qdd(6),qdd(7));
                
   robot.Jdqd_SO_q = Jdqd_b_ft_SO_q(q(1),q(2),q(3),q(4),q(5),q(6),q(7),...
                  qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7)); 
             
   robot.Jdqd_SO_v = Jdqd_b_ft_SO_v(q(1),q(2),q(3),q(4),q(5),q(6),q(7),...
                 qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7));  

   robot.Jdqd_SO_vq = Jdqd_b_ft_SO_vq(q(1),q(2),q(3),q(4),q(5),q(6),q(7),...
                  qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7));  
          
   robot.J_FO_q = J_back_ft_FO_q(q(1),q(2),q(3),q(4),q(5),q(6),q(7));
   robot.JT_FO_q = JT_back_ft_FO_q(q(1),q(2),q(3),q(4),q(5),q(6),q(7));          
             
else
    disp('Wrong Foot Selected!')
end

end

