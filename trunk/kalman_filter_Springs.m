
function [] = kalman_filter(movfilestr)

%parameters
steps = 100; % [50] number of frames to track, [-1] => all
startframe = 1; % [380] [30] starting frame
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
   steps = size(mov,2)-startframe+1;
end

img = mov(framenum).cdata;

%Define the Color Arrays
cArray = ['bo', 'g+', 'r*', 'ys'; 
            'yo', 'c+', 'm*', 'ys';
            'yo', 'c+', 'b*', 'ys';
            'yo', 'c+', 'c*', 'ys'];

psets(numsets) = spring_particle_set;
for i = 1:numsets
   
   psets(i) = spring_particle_set(steps);
   
   %Set Spring Coefficient and Spring distance threshold
   psets(i).springK = 0.1;
   psets(i).sT = 70;
   
   psets(i) = psets(i).initialize(img, cArray(i,:));
end


for i = 1:numsets
    if i > 1
        springMate = psets(i-1).target;
    else
        springMate = psets(i).target;
    end
    
    psets(i) = psets(i).set_startDistance(springMate);
end
%pause();

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
        if i > 1
            springMate = psets(i-1).target;
            psets(i) = psets(i).track_spring_target_step(s,springMate,img);
        else
            psets(i) = psets(i).track_target_step(s,img);
        end
    end
    % pause(); % interactive
    
    %Save the current image to a file
        saveas(gcf, ['ParticleFilterFrames\Frame', num2str(startframe + s - 1), '.jpg']);
        
    framenum = framenum + 1;
    
    if( s > 70 )
       disp('here') 
    end
end

%Plot the motion
close();
imshow([])
colors = ['r-', 'g-', 'b-', 'y-', 'k-'];
currentColor = 1;
for i = 1:numsets
    oldtargets = psets(i).oldtargets;
        
    plot(oldtargets(:,1),oldtargets(:,2), colors(currentColor:currentColor+1)); 
    currentColor = currentColor + 1;
end

axis([0, 700, 0, 700])


%Save the final path plot
saveas(gcf, ['ParticleFilterFrames\FinalPathNoImage.jpg']);


    
