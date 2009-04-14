
function [] = kalman_filter(movfilestr)

%parameters
steps = -1; % [50] number of frames to track, [-1] => all
startframe = 1; % [380] starting frame
framenum = startframe;

% get user input
numsets = input('\nEnter number of particle sets: ');

% load the movie
close all;
% hold on;

if steps ~= -1
   mov = aviread(movfilestr,startframe:(startframe+steps-1));
else
   mov = aviread(movfilestr);
   steps = size(mov,2);
end

img = mov(framenum).cdata;

psets(numsets) = particle_set;
for i = 1:numsets
   psets(i) = particle_set(steps);
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
    % pause(); % interactive
    
    framenum = framenum + 1;
end


    
