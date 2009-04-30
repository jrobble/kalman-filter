function [] = test_sobel(movfilestr)

mov = aviread(movfilestr);
frame = mov(100);
img = frame.cdata;
% img = imread('barcode.jpg');
close all;

gimg = rgb2gray(img);
horz = fspecial('sobel');   % create a 3-by-3 horizontal Sobel filter
vert = horz';               % create a 3-by-3 horizontal Sobel filter
Gx = filter2(horz,gimg);     % carries out filtering using the Sobel filter
Gy = filter2(vert,gimg);     % carries out filtering using the Sobel filter

% figure; imshow(mat2gray(Gx));
% figure; imshow(mat2gray(Gy));

imgheight = size(img,1);
imgwidth = size(img,2);

Sthresh = 100;
gx = zeros(imgheight,imgwidth); % range [-1,+1]
gy = zeros(imgheight,imgwidth); % range [-1,+1]
for y = 1:imgheight
   for x = 1:imgwidth
       S = sqrt(Gx(y,x)^2 + Gy(y,x)^2);
       if S > Sthresh || S == 0
           gx(y,x) = 0;
           gy(y,x) = 0;
       else
           gx(y,x) = Gx(y,x)/S;
           gy(y,x) = Gy(y,x)/S;
       end
   end
end

numbins = 8; % per component
[gdist,gcenters] = hist([-1:1],numbins);

% determine bins and create K channels
% ghist = zeros(imgheight,imgwidth);
gxchans = zeros(imgheight,imgwidth,numbins);
gychans = zeros(imgheight,imgwidth,numbins);
for y = 1:imgheight
   for x = 1:imgwidth
       [gxdiff,gxbin] = min(abs(gcenters-gx(y,x)));
       [gydiff,gybin] = min(abs(gcenters-gy(y,x)));
       % gbin = numbins*(gybin-1) + (gxbin-1) + 1; % range [1:numbins^2]
       % ghist(y,x) = gbin;
       gxchans(y,x,gxbin) = 1;
       gychans(y,x,gybin) = 1;
   end
end

igxchans = zeros(imgheight,imgwidth,numbins);
igychans = zeros(imgheight,imgwidth,numbins);
% calculate integral histogram distribution for channel
for x = 1:imgwidth
    gxprevs = zeros(1,numbins);
    gyprevs = zeros(1,numbins);
    for y = 1:imgheight
        gxs(1:numbins) = gxprevs(1:numbins) + reshape(gxchans(y,x,1:numbins),1,numbins);
        gys(1:numbins) = gyprevs(1:numbins) + reshape(gychans(y,x,1:numbins),1,numbins);
        if x-1 == 0
            gxii = gxs;
            gyii = gys;
        else
            gxii(1:numbins) = reshape(igxchans(y,x-1,1:numbins),1,numbins) + gxs(1:numbins); 
            gyii(1:numbins) = reshape(igychans(y,x-1,1:numbins),1,numbins) + gys(1:numbins); 
        end
        igxchans(y,x,1:numbins) = gxii(1:numbins);
        igychans(y,x,1:numbins) = gyii(1:numbins);
        gxprevs = gxs;
        gyprevs = gys;
    end
end

% DEBUG
% igxchans(1:10,1:10,1)
% igxchans(1:10,1:10,2)
% igxchans(1:10,1:10,3)
% igxchans(1:10,1:10,4)

% get rectangular region
% wait until user double-clicks region to resume
figure;
imshow(img);
h = imrect(gca);
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(h,fcn);
pos = wait(h); % [x1,y1,width,height], [x1,y1] denotes upper left point
delete(h);
rectangle('Position',pos,'LineWidth',2,'EdgeColor','r');

gxmodel = calc_quad_feature_model(pos,igxchans,numbins);
% gxmodel
gymodel = calc_quad_feature_model(pos,igychans,numbins);
% gymodel
tk = [gxmodel,gymodel];

while(true)
    % get rectangular region
    % wait until user double-clicks region to resume
    h = imrect(gca,pos);
    fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
    setPositionConstraintFcn(h,fcn);
    setResizable(h,false);
    pos = wait(h); % [x1,y1,width,height], [x1,y1] denotes upper left point
    delete(h);

    gxmodel = calc_quad_feature_model(pos,igxchans,numbins);
    % gxmodel
    gymodel = calc_quad_feature_model(pos,igychans,numbins);
    % gymodel
    hk = [gxmodel,gymodel];
    
    sd = 0.5;
    rho = (tk-hk).^2;
    rho = sum(rho,2);
    rho = sum(rho);
    rho = sqrt(rho);
    p = exp(-(rho.^2)/sd^2);
    p  
end


function [fmodel] = calc_quad_feature_model(pos,iimg,numf)
imgheight = size(iimg,1);
imgwidth = size(iimg,2);

% divide into equal-sized quadrants
width  = round(pos(3)/2);
height = round(pos(4)/2);
pos = round(pos);

q1pos = [pos(1)+width,   pos(2)+height,   width, height];
q2pos = [pos(1)+2*width, pos(2)+height,   width, height];
q3pos = [pos(1)+width,   pos(2)+2*height, width, height];
q4pos = [pos(1)+2*width, pos(2)+2*height, width, height];
qpos = [q1pos;q2pos;q3pos;q4pos];

% check region validity
for i = 1:size(qpos,2)
    if qpos(i,1) < 2 || qpos(i,2) < 2 || ...
       qpos(i,1) - width > imgwidth - 1 || qpos(i,2) - height > imgheight - 1 
        fmodel = [];
        return;
    end
end

% update region bounds if necessary
for i = 1:size(qpos,1)
    if qpos(i,1) > imgwidth
       qpos(i,2) = qpos(1,2) - (qpos(i,1) - imgwidth); 
       qpos(i,1) = imgwidth;
    end
    if qpos(i,2) > imgheight
       qpos(i,4) = qpos(1,4) - (qpos(i,2) - imgheight); 
       qpos(i,2) = imgheight; 
    end
end

q1pos(1:4) = qpos(1,1:4);
q2pos(1:4) = qpos(2,1:4);
q3pos(1:4) = qpos(3,1:4);
q4pos(1:4) = qpos(4,1:4);

q1sum = iimg(q1pos(2)-1,q1pos(1)-1,:);
q2sum = iimg(q2pos(2)-1,q2pos(1)-1,:);
q3sum = iimg(q3pos(2)-1,q3pos(1)-1,:);
q4sum = iimg(q4pos(2)-1,q4pos(1)-1,:);

q1k = q1sum;
q2k = q2sum-q1sum;
q3k = q3sum-q1sum;
q4k = (q4sum+q1sum)-(q2sum+q3sum); % contains no excess

% get rid of upper left data outside of ROI
if pos(1) > 1 % not along left side
    q1k(1:numf) = q1k(1:numf) - iimg(q1pos(2)-1,        q1pos(1)-q1pos(3), 1:numf);
    q3k(1:numf) = q3k(1:numf) - iimg(q3pos(2)-1,        q3pos(1)-q3pos(3), 1:numf);
    q3k(1:numf) = q3k(1:numf) + iimg(q3pos(2)-q3pos(4), q3pos(1)-q3pos(3), 1:numf);
end
if pos(2) > 1 % not along upper side
    q1k(1:numf) = q1k(1:numf) - iimg(q1pos(2)-q1pos(4), q1pos(1)-1,        1:numf);
    q2k(1:numf) = q2k(1:numf) - iimg(q2pos(2)-q2pos(4), q2pos(1)-1,        1:numf);
    q2k(1:numf) = q2k(1:numf) + iimg(q2pos(2)-q2pos(4), q2pos(1)-q2pos(3), 1:numf);
end
if pos(1) > 1 && pos(2) > 1
    q1k(1:numf) = q1k(1:numf) + iimg(q1pos(2)-q1pos(4), q1pos(1)-q1pos(3), 1:numf);
end

tmpq1k(1:numf) = q1k(:,:,1:numf)/(q1pos(3)*q1pos(4));
tmpq2k(1:numf) = q2k(:,:,1:numf)/(q2pos(3)*q2pos(4));
tmpq3k(1:numf) = q3k(:,:,1:numf)/(q3pos(3)*q3pos(4));
tmpq4k(1:numf) = q4k(:,:,1:numf)/(q4pos(3)*q4pos(4));

q1k = tmpq1k;
q2k = tmpq2k;
q3k = tmpq3k;
q4k = tmpq4k;

fmodel = [q1k;q2k;q3k;q4k];
