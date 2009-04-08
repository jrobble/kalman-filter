function [roimean,roistd] = circ_roi(movfilestr)

% image setup
[mov] = aviread(movfilestr);
img = mov(1).cdata;
imagesc(img);
axis square;

% roisetup
t = 0:pi/20:2*pi;
R0 = 5; x0 = 50; y0 = 50;
xi = R0*cos(t)+x0;
yi = R0*sin(t)+y0;
LineHandler = line(xi,yi,'LineWidth',3,'Color',[.8 0 0]);

% calc. roi stat.
roimask = poly2mask(xi,yi, size(img,1),size(img,2));
pr_r = find(roimask);
img(pr_r)
roimean = mean(img(pr_r));



end

