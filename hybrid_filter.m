
function [] = hybrid_filter(movfilestr)

%parameters
steps = -1; % [50] number of frames to track, [-1] => all
startframe = 100; % [125][200][690][60] starting frame
framenum = startframe;
numedgebins = 8; % per edge component

% get user input
numsets = input('\nEnter number of particle sets: ');

close all;

if steps ~= -1
   mov = aviread(movfilestr,startframe:(startframe+steps-1));
else
   mov = aviread(movfilestr);
   steps = size(mov,2)-startframe+1;
end

img = mov(framenum).cdata;
imshow(img);

iimg = create_color_integral(img);
[igxchans,igychans] = create_edge_integrals(img,numedgebins);
psets(numsets) = hybrid_particle_set;
for i = 1:numsets
   psets(i) = hybrid_particle_set(steps);
   psets(i) = psets(i).initialize(img,iimg,numedgebins,igxchans,igychans);
end

pause();

hold on;
imgh = image(img);

s = 0; islost = false;
while s <= steps && ~islost
    fprintf('Step %d\n',s); % DEBUG
    
    % show current image
    img = mov(framenum).cdata;
    set(imgh,'CData',img); 
    drawnow;
    
    % update all particle sets
    iimg = create_color_integral(img);
    [igxchans,igychans] = create_edge_integrals(img,numedgebins);
    for i = 1:numsets
        [psets(i),islost] = psets(i).track_target_step(s,img,iimg,igxchans,igychans);
    end
    % pause(); % interactive
    
    framenum = framenum + 1;
    s = s + 1;
end


% create integral image
function [iimg] = create_color_integral(img)
imgheight = size(img,1);
imgwidth = size(img,2);
dimg = double(img);
iimg = dimg;

for x = 1:imgwidth
    prevs = [0,0,0];
    for y = 1:imgheight
        s(1) = prevs(1) + dimg(y,x,1);
        s(2) = prevs(2) + dimg(y,x,2);
        s(3) = prevs(3) + dimg(y,x,3);
        if x-1 == 0
            ii = s;
        else
            ii(1) = iimg(y,x-1,1);
            ii(2) = iimg(y,x-1,2);
            ii(3) = iimg(y,x-1,3);
            ii = ii + s; 
        end
        iimg(y,x,1:3) = ii(1:3);
        prevs = s;
    end
end


% create edge orientation channel integral images
function [igxchans,igychans] = create_edge_integrals(img,numbins)
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
    
