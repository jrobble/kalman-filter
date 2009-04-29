function [img] = test_color_rect(movfilestr)

%parameters
steps = -1; % [50] number of frames to track, [-1] => all
startframe = 125; % [380] [30] starting frame
framenum = startframe;

close all;

if steps ~= -1
   mov = aviread(movfilestr,startframe:(startframe+steps-1));
else
   mov = aviread(movfilestr);
   steps = size(mov,2)-startframe+1;
end

img = mov(framenum).cdata;

% create integral image
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

figure; imshow(img); % DEBUG
dimg(1:10,1:10,1)
iimg(1:10,1:10,1)
dimg(1:10,1:10,2)
iimg(1:10,1:10,2)
dimg(1:10,1:10,3)
iimg(1:10,1:10,3)


%{
% get rectangular region
% wait until user double-clicks region to resume
figure;
imshow(img);
h = imrect(gca);
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(h,fcn);
pos = wait(h); % [x1,y1,width,height], [x1,y1] denotes upper left point
delete(h);

% divide into equal-sized quadrants
width  = round(pos(3)/2);
height = round(pos(4)/2);
pos = round(pos);

q1pos = [pos(1),pos(2),width,height];
q2pos = [pos(1)+width,pos(2),width,height];
q3pos = [pos(1),pos(2)+height,width,height];
q4pos = [pos(1)+width,pos(2)+height,width,height];

% c = rectangle('Position',pos,'LineWidth',2,'EdgeColor','c');
q1r = rectangle('Position',q1pos,'LineWidth',1,'EdgeColor','r');
q2r = rectangle('Position',q2pos,'LineWidth',1,'EdgeColor','g');
q3r = rectangle('Position',q3pos,'LineWidth',1,'EdgeColor','b');
q4r = rectangle('Position',q4pos,'LineWidth',1,'EdgeColor','y');

q1k = calc_color_model_region(q1pos,img);
q2k = calc_color_model_region(q2pos,img);
q3k = calc_color_model_region(q3pos,img);
q4k = calc_color_model_region(q4pos,img);
tk = [q1k;q2k;q3k;q4k];

% DEBUG
figure;
imshow(img);
% c = rectangle('Position',pos,'LineWidth',1,'EdgeColor','c');
q1r = rectangle('Position',q1pos,'EdgeColor','none','FaceColor',q1k/255);
q2r = rectangle('Position',q2pos,'EdgeColor','none','FaceColor',q2k/255);
q3r = rectangle('Position',q3pos,'EdgeColor','none','FaceColor',q3k/255);
q4r = rectangle('Position',q4pos,'EdgeColor','none','FaceColor',q4k/255);

while(true)
    % select another rectangle region of same size for comparison
    h = imrect(gca,pos);
    fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
    setPositionConstraintFcn(h,fcn);
    setResizable(h,false);
    hpos = wait(h);
    delete(h);

    % sum color componenet diffs for each region
    hk = calc_quad_color_model(hpos,img);
    sig = 10;
    
    % rho = sqrt(sum((tk-hk).^2,2));

    rho = (tk-hk).^2;
    rho = sum(rho,2);
    rho = sum(rho);
    rho = sqrt(rho);
    
    p = exp(-(rho.^2)/sig^2);
    
    % DEBUG
    tk
    hk
    rho
    p
end



function [cmodel] = calc_quad_color_model(pos,img)
% divide into equal-sized quadrants
width  = round(pos(3)/2);
height = round(pos(4)/2);
pos = round(pos);

q1pos = [pos(1),pos(2),width,height];
q2pos = [pos(1)+width,pos(2),width,height];
q3pos = [pos(1),pos(2)+height,width,height];
q4pos = [pos(1)+width,pos(2)+height,width,height];

q1k = calc_color_model_region(q1pos,img);
q2k = calc_color_model_region(q2pos,img);
q3k = calc_color_model_region(q3pos,img);
q4k = calc_color_model_region(q4pos,img);
cmodel = [q1k;q2k;q3k;q4k];


function [rcmodel] = calc_color_model_region(rpos,img)
area = rpos(3)*rpos(4);
rcmodel = zeros(1,3);
img = double(img);

for y = rpos(2):rpos(2)+rpos(4)-1
    for x = rpos(1):rpos(1)+rpos(3)-1
      rcmodel(1) = rcmodel(1) + img(y,x,1);
      rcmodel(2) = rcmodel(2) + img(y,x,2);
      rcmodel(3) = rcmodel(3) + img(y,x,3);
   end
end

rcmodel = round(rcmodel / area);
%}