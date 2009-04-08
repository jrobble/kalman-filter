function [] = test_students()

% mplay('Students.avi'); % works
% [mov] = aviread('Students.avi', [1:200]);
% img = mov(100).cdata;
% imshow(img);

% aviinfo('Students.avi')
% [video, audio] = mmread('Students.avi',[1:1461]);
% mplay(video.frames);

play_matlab_fbf();

% split();


function [] = play_mmread_fbf() 

aviinfo('Students_1.avi')
[video, audio] = mmread('Students_1.avi');
for i = 1:size(video.frames,2)
   imshow(video.frames(i).cdata);
   drawnow;
end


function [] = play_matlab_fbf() 

aviinfo('Students_1.avi')
[mov] = aviread('Students_1.avi');
size(mov,2)
for i = 1:size(mov,2)
   imshow(mov(i).cdata);
   drawnow;
end


function [] = split()
    % alter path
    path(path,'mmread')

    splitframe = 1461; % first frame of second half
    numframes = 3214;

    % use same compression and frame rate as original video
    newmov = avifile('Students_1.avi', ...
        'compression', 'Cinepak','fps',29.9700,'quality',100);
    [video, audio] = mmread('Students.avi',1:(splitframe-1));

    % part 1
    for i = 1:(splitframe-1)
       frame = video.frames(i);
       newmov = addframe(newmov,frame);
    end

    % part 2
    newmov = avifile('Students_2.avi', ...
        'compression', 'Cinepak','fps',29.9700,'quality',100);
    [video, audio] = mmread('Students.avi',splitframe:numframes);
    
    for i = 1:(numframes-splitframe)
       frame = video.frames(i);
       newmov = addframe(newmov,frame);
    end

    clear all;