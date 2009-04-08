% RGB space is fine since we're not worried about lighting conditions
% (otherwise use HSV space).
% Update rule uses Brownian particle motion.
classdef particle_set
    
    properties
        % parameters
        numbins = 3; % [3] number of histogram bins to use for each color component
        usehsv = 1; % [1]

        rad = 5;  % [5] radius of circular ROI
        var = 0.2; % [0.02] variance of Gaussian used to calculate probabilities (lower is more selective)
        scalex = 15; % approx. max Brownian particle x movement (pos/neg)
        scaley = 15; % approx. max Brownian particle y movement (pos/neg)

        winrad = 20; % radius of small circular window around best particle in "robust mean" method  
        numpts = 20; % [200] number of particles in system

        % USC helicopter video:
        % frame 350 - helicopter low in tree before occlusion by man walking in front
        % frame 200 - basic floating helicopter
        % numbins = 3; % [8][3]
        % usehsv = 0;
        % rad = 5;  % [10][5]
        % var = 0.01; % [0.001][0.01]
        % scalex = 15; % [15][5]
        % scaley = 15; % [15][5]
        % winrad = 20; % [20]
        % numpts = 200; 

        % ACV student video:
        % frame 200 - Brian standing

        % tracking variables
        graphics = []; % set of graphic handles
        steps = []; % number of steps to track
        sk = [];
        oldsk = [];
        oldtargets = [];
        target = [];
        q = [];
    end
    
    methods
        
        % constructor must be able to handle the no input case for object array pre-allocation
        function obj = particle_set(steps)
            if nargin > 0
                obj.steps = steps;
            end
        end
        
        function obj = initialize(obj,img)
            % choose region of interest
            % use defaults, let user select circle center and radius
            % [circhandle,circcenter,circrad]=draw_circle([],[],[],[],img);

            % create initial set of particles
            %{
            imshow(img);
            fprintf(1,'Select target point.\n\n');
            [x,y] = ginput(1);
            target = [x,y];

            fprintf(1,'Select points in particle population, then press Enter.\n');
            [x,y] = ginput();
            numpts = size(x,1);
            %}

            % select uniformly distributed random points w/in circle around target point
            fprintf(1,'Select ROI. Target point is the center.\n\n');
            [circhandle,circcenter,circrad]=draw_circle([],[],[],[],img);
            % obj.rad = circrad; % DEBUG
            obj.target = circcenter;
            x = zeros(1,obj.numpts);
            y = zeros(1,obj.numpts);
            count = 1;
            while count <= obj.numpts
               tmpx = unifrnd(circcenter(1)-circrad,circcenter(1)+circrad);
               tmpy = unifrnd(circcenter(2)-circrad,circcenter(2)+circrad);
               tmpdist = sqrt((tmpx-circcenter(1))^2 + (tmpy-circcenter(2))^2);
               if tmpdist <= circrad
                  % fprintf(1,'+ tmpx: %d tmpy: %d\n',tmpx,tmpy); % DEBUG
                  x(count) = tmpx;
                  y(count) = tmpy;
                  count = count + 1;
               end
            end

            % sk = [x,y,b,c]; 
            % [xcoord, ycoord, 
            %  probability, cumulative probability, 
            %  Brownian x, Brownian y, Brownian dx, Brownian dy]
            obj.sk = zeros(obj.numpts,6);
            obj.sk(:,1) = x(:);
            obj.sk(:,2) = y(:);
            obj.sk(:,5) = 1;
            obj.sk(:,6) = 1;
            obj.sk(:,7) = 1;
            obj.sk(:,8) = 1;

            fprintf(1,'Selected points:\n');
            for i = 1:obj.numpts
               fprintf(1,'[%d,%d]\n',obj.sk(i,1),obj.sk(i,2)); 
            end
            fprintf(1,'\n');

            % calculate target color distribution
            % fprintf(1,'Calculating target color distribution, [%d,%d] ...\n',target(1),target(2)); % DEBUG
            if obj.usehsv == 1
                img = rgb2hsv(img);
            end
            obj.q = obj.get_color_distribution(obj.target,img);

            % step through alg.
            obj.oldtargets = zeros(obj.steps,2);
        end
        
        
        function obj = track_target_step(obj,step,img)
            % clear previous graphics
            for i = 1:size(obj.graphics,2)
               delete(obj.graphics(i)); 
            end
            
            % mark points (particles)
            gcount = 1;
            for i = 1:size(obj.oldsk,1)
                obj.graphics(gcount) = plot(obj.oldsk(:,1),obj.oldsk(:,2),'b+');
                gcount = gcount + 1;
            end
            obj.graphics(gcount) = plot(obj.sk(:,1),obj.sk(:,2),'g+');
            gcount = gcount + 1;

            % mark targets (tracking)
            obj.oldtargets(step,:) = obj.target(:);
            for i = 1:step-1
               obj.graphics(gcount) = plot(obj.oldtargets(i,1),obj.oldtargets(i,2),'r*');
               gcount = gcount + 1;
            end

            obj.graphics(gcount) = plot(obj.target(1),obj.target(2),'ys'); 
            % drawnow; % bottleneck 

            if obj.usehsv == 1
                img = rgb2hsv(img);
            end

            % calculate particle probabilities
            for i = 1:obj.numpts
                % fprintf(1,'Calculating particle color distribution, [%d,%d] ...\n',sk(i,1),sk(i,2)); % DEBUG
                [py] = obj.get_color_distribution(obj.sk(i,:),img);

                % calculate Bhattacharyya coefficient (ro)
                ro = sum(sqrt(py .* obj.q));

                % calculate probability
                const = 1/(sqrt(2*pi*obj.var));
                pow = -(1-ro)/(2*obj.var);
                b = const*exp(pow);
                obj.sk(i,3) = b;

                % calculate cumulative probability
                c = sum(obj.sk(1:i,3));
                obj.sk(i,4) = c;
                
                % DEBUG
                h = plot(obj.sk(i,1),obj.sk(i,2),'w+');
                obj.sk(i,3)
                pause();
                delete(h);
            end


            % normalize the probabilities
            normfactor = sum(obj.sk(:,3));
            obj.sk(:,3) = obj.sk(:,3)/normfactor;
            obj.sk(:,4) = obj.sk(:,4)/normfactor;

            % DEBUG
            %{
            for i = 1:numpts
               fprintf(1,'sk(%d) b: %d, c: %d\n',i,sk(i,3),sk(i,4));
            end
            %}

            % determine the target pt which defines the next object position
            % use 1 of 3 approaches 

            % 1. "weighted mean" method 
            % target = [sum(sk(:,3) .* sk(:,1)), ...
            %           sum(sk(:,3) .* sk(:,2))];
            % target % DEBUG

            % 2. "best particle" method ...
            [val,index] = max(obj.sk(:,3));
            best(1:8) = obj.sk(index,1:8);
            obj.target(1:2) = best(1:2);

            % 3. "robust mean" method (weighted mean in small circular window around best particle)
            % plot(best(1),best(2),'w+'); % DEBUG
            %{
            obj.target = zeros(1,2);
            normfactor = 0;
            count = 0;
            for i = 1:obj.numpts
                if (obj.sk(i,1)-best(1))^2+(obj.sk(i,2)-best(2))^2 < obj.winrad^2
                    normfactor = normfactor + obj.sk(i,3);
                    obj.target(1) = obj.target(1) + obj.sk(i,3) * obj.sk(i,1);
                    obj.target(2) = obj.target(2) + obj.sk(i,3) * obj.sk(i,2);
                    count = count + 1;
                end
            end
            if count == 0 % should never happen - the best particle should always be in the window
               error('Robust mean failed: no neighbors around best particle.'); 
            end
            obj.target(:) = obj.target(:) / normfactor;
            % obj.target % DEBUG
            %}


            % resample particles (can choose repeat particles)
            % fprintf(1,'\nResample particles ...\n'); % DEBUG
            sk0 = zeros(obj.numpts,size(obj.sk,2));
            for i = 1:obj.numpts
                r = rand(); % random num normally dist. in interval [0,1]
                % find smallest particle cumulative density >= r
                for j = 1:obj.numpts
                   if obj.sk(j,4) >= r
                      sk0(i,1:4) = obj.sk(j,1:4); % particle inherits new location
                      sk0(i,5:8) = obj.sk(i,5:8); % particle maintains Brownian paramters (based on index) 
                      break;
                   end
                end
                if ~any(sk0(i))
                   obj.sk
                   error('No particle chosen for r: %d\nPossible loss of target.',r); 
                end
                % fprintf(1,'r: %d, sk0(%d): [%d,%d,%d,%d]\n', ...
                %    r,i,sk0(i,1),sk0(i,2),sk0(i,3),sk0(i,4)); % DEBUG
            end
            % sk0 % DEBUG


            % DEBUG - mark points (particles) before update rule
            xmin = min(sk0(:,1));
            xmax = max(sk0(:,1));
            ymin = min(sk0(:,2));
            ymax = max(sk0(:,2));
            box = [[xmin,ymin]; ...
                   [xmax,ymin]; ...
                   [xmax,ymax]; ...
                   [xmin,ymax]; ...
                   [xmin,ymin]];
            %{
            imshow(img,'InitialMagnification','fit');
            for i = 1:size(sk0,1)
                plot(sk0(:,1),sk0(:,2),'w+'); 
            end
            plot(box(:,1),box(:,2),'y');
            pause();
            %}

            % update Brownian motion scale factors
            % obj.scalex = (xmax - xmin)/2 + 2;
            % obj.scaley = (ymax - ymin)/2 + 2;
            
            maxb = best(3);
            maxb

            % DEBUG 
            %{
            x = randi(100);
            y = randi(100);
            for i = 1:numpts
                sk0(i,:) = [x,y,0,0,sk(i,5),sk(i,6),sk(i,7),sk(i,8)];
            end
            plot(x,y,'w+'); 
            pause();
            %}

            % spread the states according to the Brownian update rule
            sk1 = zeros(obj.numpts,size(obj.sk,2));
            for i = 1:obj.numpts
                sk1(i,5) = exp((-1/4) * (sk0(i,5) + 1.5*sk0(i,7)));
                sk1(i,6) = exp((-1/4) * (sk0(i,6) + 1.5*sk0(i,8)));
                sk1(i,7) = exp((-1/4) * sk0(i,7));
                sk1(i,8) = exp((-1/4) * sk0(i,8));
                % add the effect of a multivariate Gaussian random variable
                sk1(i,5) = sk1(i,5) + rand();
                sk1(i,6) = sk1(i,6) + rand();
                sk1(i,7) = sk1(i,7) + rand();
                sk1(i,8) = sk1(i,8) + rand();
                % put in proper range
                sk1(i,1) = (sk1(i,5)-sk0(i,5)) * obj.scalex + sk0(i,1);
                sk1(i,2) = (sk1(i,6)-sk0(i,6)) * obj.scaley + sk0(i,2);
            end

            % sk1 % DEBUG
            obj.oldsk = obj.sk;
            obj.sk = sk1;
        end


        % get color histogram distribution for a pixel
        % p - pixel
        % rad - ROI circle radius
        % img - image
        function [py] = get_color_distribution(obj,p,img)
            % calculate color histogram profile
            imgheight = size(img,1);
            imgwidth = size(img,2);
            centerx = p(1);
            centery = p(2);

            % determine boundaries of window to check
            ymin = max(1,floor(centery-obj.rad));
            ymax = min(imgheight,ceil(centery+obj.rad));
            xmin = max(1,floor(centerx-obj.rad));
            xmax = min(imgwidth,ceil(centerx+obj.rad));

            % subimg = uint8(zeros(imgheight,imgwidth,3)); % DEBUG
            area = floor(pi*obj.rad^2);
            innerpix = zeros(area,3);
            count = 1;
            for y = ymin:ymax
                for x = xmin:xmax
                    % determine if pixel within circle ROI
                    if (x-centerx)^2+(y-centery)^2 <= obj.rad^2
                        % subimg(y,x,:) = img(y,x,:); % DEBUG
                        innerpix(count,:) = img(y,x,:);
                        count = count + 1;
                    end

                end
            end
            count = count - 1;
            innerpix = innerpix(1:count,:);
            % fprintf(1,'innerpix: %d\n',count); % DEBUG

            % figure(2); % DEBUG
            % imshow(subimg); %DEBUG

            % calculate color histogram centers for ROI
            [rdist,rcenters] = hist(innerpix(:,1),obj.numbins);
            [gdist,gcenters] = hist(innerpix(:,2),obj.numbins);
            [bdist,bcenters] = hist(innerpix(:,3),obj.numbins);
            ccenters = [rcenters;gcenters;bcenters];
            % ccenters % DEBUG

            % calculate color distribution at location y w/in the ROI 
            % for all m = 8x8x8 = 512 bins (if numbins == 8)
            py = zeros(1,obj.numbins^3);
            normfactor = 0;

            for y = ymin:ymax
                for x = xmin:xmax
                    % determine if pixel within circle ROI
                    if (x-centerx)^2+(y-centery)^2 <= obj.rad^2
                        % absolute Euclidean distance divided by radius
                        e = pdist([[x,y];[centerx,centery]]);
                        e = abs(e)/obj.rad;

                        k = 0;
                        if e < 1 % should always be true
                            k = 1-e^2;  
                        else
                            error('e >= 1, %d >= 1',e);
                        end

                        % get color histogram profile
                        p(1:3) = img(y,x,:);
                        hbins = obj.get_color_profile(p,ccenters);

                        % determine bin index (like applying Kronecker delta)
                        u = (hbins(3) * obj.numbins^2) + (hbins(2) * obj.numbins) + hbins(1) + 1;
                        % hbins % DEBUG
                        % u % DEBUG

                        % calculate the norm factor once
                        normfactor = normfactor + k;

                        % calculate the probability density
                        py(u) = py(u) + k;

                        % error('DONE'); % DEBUG
                    end

                end
            end
            % fprintf(1,'>> ubin [%d,%d,%d]: %d, tcount: %d\\%d\n', ...
            %    rbin,gbin,bbin,ucount,tcount,count); % DEBUG

            % apply norm factor to distribution
            normfactor = 1/normfactor;
            py = normfactor * py;
            % fprintf(1,'>> sum(py): %d\n',sum(py)); % DEBUG
        end


        % get color histogram profile for a pixel
        % p - pixel
        % ccenters - color histogram centers
        function [cprofile] = get_color_profile(obj,p,ccenters)
            % ccenters % DEBUG
            % p % DEBUG
            p = double(p);
            [mindiffr,posr] = min(abs(p(1)-ccenters(1,:)));
            [mindiffg,posg] = min(abs(p(2)-ccenters(2,:)));
            [mindiffb,posb] = min(abs(p(3)-ccenters(3,:)));
            cprofile = [posr-1,posg-1,posb-1];
            % cprofile % DEBUG
        end
        
    end % methods
end % class

