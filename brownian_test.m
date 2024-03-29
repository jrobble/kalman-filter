function [] = brownian_test()

% scale and center
%{
scaley = 5;
scalex = 10;
centerx = 50;
centery = 100;
%}
scaley = 1;
scalex = 1;
centerx = 1;
centery = 1;

steps = 100;
sk = zeros(steps,4); % x, y, dx, dy
sk(1,:) = [1,1,1,1]; % tends to center around (1,1) in range [[0,2],[0,2]]
skreal = zeros(steps,2);
skreal(1,:) = [centerx,centery];

for i = 1:(steps-1)
    sk(i+1,1) = exp((-1/4) * (sk(i,1) + 1.5*sk(i,3)));
    sk(i+1,2) = exp((-1/4) * (sk(i,2) + 1.5*sk(i,4)));
    sk(i+1,3) = exp((-1/4) * sk(i,3));
    sk(i+1,4) = exp((-1/4) * sk(i,4));
    
    % add the effect of a multivariate Gaussian random variable
    % (rand()*(maxx-minx)+minx);
    sk(i+1,1) = sk(i+1,1) + rand();
    sk(i+1,2) = sk(i+1,2) + rand();
    sk(i+1,4) = sk(i+1,3) + rand();
    sk(i+1,5) = sk(i+1,4) + rand();
    
    fprintf(1,'step: %d, sk[%d,%d,%d,%d]\n',i,sk(i,1),sk(i,2),sk(i,3),sk(i,4));
    
    % put in proper range
    skreal(i+1,1) = (sk(i+1,1)-sk(i,1)) * scalex + centerx;
    skreal(i+1,2) = (sk(i+1,2)-sk(i,2)) * scaley + centery;
    
    fprintf(1,'step: %d, skreal[%d,%d]\n',i,skreal(i,1),skreal(i,2));
end
        
% figure(1);
% plot([1:steps],sk(:,1),'r');
   
figure(2);
plot(skreal(:,1),skreal(:,2),'b');
% axis([0,2,0,2]);
xlabel('x');
ylabel('y');

figure(3);
plot3([1:steps],skreal(:,1),skreal(:,2),'b');
% axis([0,100,0,2,0,2]);
grid;
xlabel('steps');
ylabel('x');
zlabel('y');

figure(4);
hist(skreal(:,1));
% xlim([0,2]);
title('x histogram');

figure(5);
hist(skreal(:,2));
% xlim([0,2]);
title('y histogram');

