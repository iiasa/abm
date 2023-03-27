function [epsilon_Y_EA,epsilon_E,epsilon_I] = epsilon(C)
L=chol(C,'lower');
epsilon=randn(3,1);
epsilon_Y_EA=L(1,:)*epsilon;
epsilon_E=L(2,:)*epsilon;
epsilon_I=L(3,:)*epsilon;
end

