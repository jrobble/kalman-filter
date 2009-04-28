% plot a circle
function [handle] = plot_circle(rad,xcenter,ycenter) 
    t = linspace(0,2*pi,1000);
    x = rad * cos(t) + xcenter;                  
    y = rad * sin(t) + ycenter;
    handle = plot(x,y,'w-');
end