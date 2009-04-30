
function [] = hybrid_filter(movfilestr)

%parameters
steps = -1; % [50] number of frames to track, [-1] => all
startframe = 1; % [125][200][690] starting frame
framenum = startframe;

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

iimg = create_integral_img(img);
psets(numsets) = hybrid_particle_set;
for i = 1:numsets
   psets(i) = hybrid_particle_set(steps);
   psets(i) = psets(i).initialize(img,iimg);
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
    iimg = create_integral_img(img);
    for i = 1:numsets
        [psets(i),islost] = psets(i).track_target_step(s,img,iimg);
    end
    % pause(); % interactive
    
    framenum = framenum + 1;
    s = s + 1;
end


% create integral image
function [iimg] = create_integral_img(img)
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

    
