function [] = test_unifrnd()

%{
hold on;
for i = 1:1000
    r = unifrnd(1,10,1,2);
    plot(r(1),r(2),'b+');
end
%}

M = [0 0];
C = [1 0; 0 1]; 
z = repmat(M,10000,1) + randn(10000,2)*C; % randn returns N by M matrix

plot(z(:,1),z(:,2),'b+');
