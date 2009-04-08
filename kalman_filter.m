
function [] = kalman_filter(movfilestr)

%parameters
steps = 50; % [50] number of frames to track
startframe = 380; % starting frame
framenum = 1;

% get user input
numsets = input('\nEnter number of particle sets: ');

% load the movie
close all;
% hold on;

mov = aviread(movfilestr,startframe:(startframe+steps));
img = mov(1).cdata;

psets(numsets) = particle_set;
for i = 1:numsets
   psets(i) = particle_set(steps);
   psets(i) = psets(i).initialize(img);
end

hold on;
imgh = image(img);
for s = 1:steps
    fprintf('Step %d\n',s); % DEBUG
    
    % show current image
    img = mov(framenum).cdata;
    set(imgh,'CData',img);
    drawnow;
    framenum = framenum + 1;
            
    % update all particle sets
    for i = 1:numsets
        psets(i) = psets(i).track_target_step(s,img);
    end
    % pause();
end


    
