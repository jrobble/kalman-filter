
function [] = hybrid_filter(movfilestr)

%parameters
steps = -1; % [50] number of frames to track, [-1] => all
startframe = 125; % [380] [30] starting frame
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

psets(numsets) = hybrid_particle_set;

imshow(img);
for i = 1:numsets
   psets(i) = hybrid_particle_set(steps);
   psets(i) = psets(i).initialize(img);
end

pause();

hold on;
imgh = image(img);
for s = 1:steps
    fprintf('Step %d\n',s); % DEBUG
    
    % show current image
    img = mov(framenum).cdata;
    set(imgh,'CData',img); 
    drawnow;
    
    % update all particle sets
    for i = 1:numsets
        psets(i) = psets(i).track_target_step(s,img);
    end
    pause(); % interactive
    
    framenum = framenum + 1;
end

    
