function [] = test()

    [mov] = aviread('Second.avi');
    img = mov(100).cdata;
    imshow(img);
    
    imshow(img);
    fprintf(1,'Select target point.\n\n');
    [x,y] = ginput(1);
    target = [x,y];
    
    fprintf(1,'Select points in particle population, then press Enter.\n');
    [x,y] = ginput();
    numpts = size(x,1);
    
    % sk = [x,y,b,c]; xcoord, ycoord, probability, cumulative probability, dx, dy
    sk = zeros(numpts,6);
    sk(:,1) = x(:);
    sk(:,2) = y(:);
    sk(:,5) = 10;
    sk(:,6) = 10;

    hold on;
    plot(target(1),target(2),'r*');
    plot(sk(:,1),sk(:,2),'g+');
    drawnow;
    
    hold on;
    plot(target(1)+10,target(2)+10,'r*');
    plot(sk(:,1)+10,sk(:,2)+10,'g+');
    drawnow;

end

