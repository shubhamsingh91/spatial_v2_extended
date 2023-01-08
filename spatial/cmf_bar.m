function out = cmf_bar(x)

out = [-skew(x(1:3)), -skew(x(4:6)) ; -skew(x(4:6)) ,zeros(3,3)];

end


