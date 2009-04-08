function [] = test_unifrnd()

hold on;
for i = 1:1000
    r = unifrnd(1,10,1,2);
    plot(r(1),r(2),'b+');
end