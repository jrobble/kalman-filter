% RGB space is fine since we're not worried about lighting conditions
% (otherwise use HSV space).
% Update rule uses Brownian particle motion.
classdef particle_set
    
    properties
        % parameters
        numbins = 3; % [3] number of histogram bins to use for each color component
        usehsv = 0; % [1]

        rad = 5;  % [5] radius of circular ROI
        var = 0.2; % [0.02] variance of Gaussian used to calculate probabilities (lower is more selective)
        origscalex = 15; % approx. max Brownian particle x movement (pos/neg)
        origscaley = 15; % approx. max Brownian particle y movement (pos/neg)

        winrad = 20; % radius of small circular window around best particle in "robust mean" method
        orignumpts = 100; % [20][200] number of particles in system
        
        % method used to determine next target pt
        % 1. "weighted mean" method 
        % 2. "best particle" method
        % 3. "robust mean" method
        targetmethod = 3;

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
        qccenters = [];
        oldnumpts = [];
        numpts = [];
        scalex = [];
        scaley = [];
        imgheight = [];
        imgwidth = [];
        islost = [];
    end
    
    methods
        
        % constructor must be able to handle the no input case for object array pre-allocation
        function obj = particle_set(steps)
            obj.islost = 1; % initially lost
            obj.numpts = obj.orignumpts;
            obj.oldnumpts = obj.numpts;
            obj.scalex = obj.origscalex;
            obj.scaley = obj.origscaley;
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
            
            obj.imgheight = size(img,1);
            obj.imgwidth = size(img,2);

            % select uniformly distributed random points w/in circle around target point
            fprintf(1,'Select ROI. Target point is the center.\n\n');
            [circhandle,circcenter,circrad]=draw_circle([],[],[],[],img);
            
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
            [obj.q,obj.qccenters] = obj.get_color_distribution(obj.target,img);
        end
        
        
        % step through alg.
        function obj = track_target_step(obj,step,img)
            if obj.usehsv == 1
                img = rgb2hsv(img);
            end

            % calculate particle probabilities
            tmpgraphics = [];
            for i = 1:obj.numpts
                % fprintf(1,'Calculating particle color distribution, [%d,%d] ...\n',sk(i,1),sk(i,2)); % DEBUG
                [py] = obj.get_color_distribution(obj.sk(i,:),img,obj.qccenters);

                if isnan(py) % DEBUG
                    obj.sk(i,:)
                    error('py is NAN');
                end
                
                % calculate Bhattacharyya coefficient (ro)
                ro = sum(sqrt(py .* obj.q));
                
                if isnan(ro) % DEBUG
                    error('ro is NAN');
                end

                % calculate probability
                const = 1/(sqrt(2*pi*obj.var));
                pow = -(1-ro)/(2*obj.var);
                b = const*exp(pow);
                obj.sk(i,3) = b;
                
                if isnan(b) % DEBUG
                    error('b is NAN');
                end

                % calculate cumulative probability
                c = sum(obj.sk(1:i,3));
                obj.sk(i,4) = c;
                
                % DEBUG
                % h = plot(obj.sk(i,1),obj.sk(i,2),'w+');
                % g = obj.plot_circle(obj.rad,obj.sk(i,1),obj.sk(i,2));
                % obj.sk(i,3)
                % pause();
                % delete(h); 
                % delete(g);
                % tmpgraphics(end+1) = h;
            end
            %pause();
            %for i = 1:size(tmpgraphics,2)
            %   delete(tmpgraphics(i)); 
            %end


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
            
            if obj.targetmethod == 1
                % 1. "weighted mean" method 
                obj.target = [sum(obj.sk(:,3) .* obj.sk(:,1)), ...
                              sum(obj.sk(:,3) .* obj.sk(:,2))];
                obj.target % DEBUG
                
            elseif obj.targetmethod == 2
                % 2. "best particle" method
                [val,index] = max(obj.sk(:,3));
                best(1:8) = obj.sk(index,1:8);
                fprintf(1,'Best particle method best: %d\n',best(3)); % DEBUG
                obj.target(1:2) = best(1:2);
                
            elseif obj.targetmethod == 3
                % 3. "robust mean" method (weighted mean in small circular window around best particle)
                % plot(best(1),best(2),'w+'); % DEBUG
                [val,index] = max(obj.sk(:,3));
                best(1:8) = obj.sk(index,1:8);
                fprintf(1,'Robust mean method best: %d\n',best(3)); % DEBUG
                
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
            end
            
            
            % DEBUG - mark points (particles) before update rule
            %{
            xmin = min(sk0(:,1));
            xmax = max(sk0(:,1));
            ymin = min(sk0(:,2));
            ymax = max(sk0(:,2));
            box = [[xmin,ymin]; ...
                   [xmax,ymin]; ...
                   [xmax,ymax]; ...
                   [xmin,ymax]; ...
                   [xmin,ymin]];
            %}
            
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
            if isnan(maxb) 
                error('maxb is NAN');
            end
            if maxb < 0.02 
                if obj.islost == 0
                    % increase particle cloud size
                    fprintf(1,'>> Possible loss of target. ');
                    fprintf(1,'=> Use larger particle cloud and best particle method.\n');
                    obj.scalex = obj.origscalex * 5;
                    obj.scaley = obj.origscaley * 5;
                    obj.targetmethod = 2; % use "best particle" method
                    obj.oldnumpts = obj.numpts;
                    % obj.numpts = obj.orignumpts;
                    obj.islost = 1;
                end
            else
                if obj.islost == 1
                    % reset particle cloud size back to default
                    fprintf(1,'>> Possible regain of target. ');
                    fprintf(1,'=> Use smaller particle cloud and robust mean method.\n');
                    obj.scalex = obj.origscalex;
                    obj.scaley = obj.origscaley;
                    obj.targetmethod = 3; % use "robust mean" method
                    
                    % might want to scan the whole image here...
                    % reset Brownian motion parameters ?
                    % increase number of particles ?
                    obj.oldnumpts = obj.numpts;
                    % obj.numpts = 2 * obj.orignumpts;
                    obj.islost = 0;
                end
            end

            
            % select resample particles (can choose repeat particles)
            % fprintf(1,'\nResample particles ...\n'); % DEBUG
            sk0 = zeros(obj.numpts,size(obj.sk,2));
            for i = 1:obj.numpts
                r = rand(); % random num normally dist. in interval [0,1]
                % find smallest particle cumulative density >= r
                found = 0; j = 1;
                while j <= obj.oldnumpts && found == 0
                   if obj.sk(j,4) >= r
                      % r
                      % h = plot(obj.sk(j,1),obj.sk(j,2),'w+'); % DEBUG
                      % pause();
                      % delete(h);
                      found = 1;
                   else
                      j = j + 1;
                   end
                end
                % if there are no good choices, choose a random particle 
                if found == 0
                   % obj.sk
                   fprintf(1,'>> Possible loss of target. ');
                   fprintf(1,'=> Arbitrary particle chosen for r: %d\n',r);
                   j = randi(obj.oldnumpts);
                end
                % set new particle parameters
                sk0(i,1:4) = obj.sk(j,1:4); % particle inherits new location
                if i <= obj.oldnumpts
                    sk0(i,5:8) = obj.sk(i,5:8); % particle maintains Brownian parameters (based on index)
                else
                    sk0(i,5:8) = obj.sk(randi(obj.oldnumpts),5:8); % if new particle, inherit random Brownian parameters
                end
                
                % fprintf(1,'r: %d, sk0(%d): [%d,%d,%d,%d]\n', ...
                %    r,i,sk0(i,1),sk0(i,2),sk0(i,3),sk0(i,4)); % DEBUG
            end
            % sk0 % DEBUG
            

            % spread the states according to the Brownian update rule
            sk1 = zeros(obj.numpts,size(obj.sk,2));
            for i = 1:obj.numpts
                isvalid = 0;
                attempts = 0;
                while isvalid == 0
                    % attempt to update
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
                    
                    % validate
                    isvalid = (sk1(i,1) > 0) && (sk1(i,1) <= obj.imgwidth) && ...
                              (sk1(i,2) > 0) && (sk1(i,2) <= obj.imgheight);
                    attempts = attempts + 1;
                    
                    % after so many attempts, just choose a random pt
                    if isvalid == 0 && attempts > 10
                        fprintf(1,'>> Too many Brownian motion parameter update attempts. ');
                        fprintf(1,'=> Choose random coordinates.\n');
                        sk1(i,1) = randi(obj.imgwidth);
                        sk1(i,2) = randi(obj.imgheight);
                        isvalid = 1;
                    end
                end
            end

            % sk1 % DEBUG
            obj.oldsk = obj.sk;
            obj.sk = sk1;
            
            % mark-up the image after determining points
            obj = obj.mark_up(step);
        end


        % get color histogram distribution for a pixel
        % p - pixel
        % rad - ROI circle radius
        % img - image
        function [py,ccenters] = get_color_distribution(obj,p,img,ccenters)
            % calculate color histogram profile
            centerx = p(1);
            centery = p(2);

            % determine boundaries of window to check
            ymin = max(1,floor(centery-obj.rad));
            ymax = min(obj.imgheight,ceil(centery+obj.rad));
            xmin = max(1,floor(centerx-obj.rad));
            xmax = min(obj.imgwidth,ceil(centerx+obj.rad));

            % calculate color histogram centers
            if nargin < 4
                fprintf(1,'>> Calculate color histogram centers.\n'); % DEBUG
                
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
            end

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
                        if e < 1 % should always be true unless random coordinate chosen
                            k = 1-e^2;  
                        else
                            % e % DEBUG
                            % error('e >= 1, %d >= 1',e);
                        end

                        % get color histogram profile
                        p(1:3) = img(y,x,:);
                        hbins = obj.get_color_profile(p,ccenters);

                        % determine bin index (like applying Kronecker delta)
                        u = (hbins(3) * obj.numbins^2) + (hbins(2) * obj.numbins) + hbins(1) + 1;

                        % calculate the norm factor once
                        normfactor = normfactor + k;

                        % calculate the probability density
                        py(u) = py(u) + k;

                        % DEBUG
                        % hold on;
                        % h = plot(x,y,'c+');
                        % u
                        % pause();
                        % delete(h);
                        
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
        
        
        % plot a circle
        function [handle] = plot_circle(obj,rad,xcenter,ycenter) 
            t = linspace(0,2*pi,1000);
            x = rad * cos(t) + xcenter;                  
            y = rad * sin(t) + ycenter;
            handle = plot(x,y,'w-');
        end
        
        
        % mark-up the image
        function obj = mark_up(obj,step)
           % clear previous graphics
            for i = 1:size(obj.graphics,2)
               delete(obj.graphics(i)); 
            end
            
            % mark points (particles)
            gcount = 1;
            obj.graphics(gcount) = plot(obj.oldsk(:,1),obj.oldsk(:,2),'b+');
            gcount = gcount + 1;
            
            obj.graphics(gcount) = plot(obj.sk(:,1),obj.sk(:,2),'g+');
            gcount = gcount + 1;

            % mark targets (tracking)
            if obj.islost == 0 % only record good targets
                obj.oldtargets(end+1,:) = obj.target(:);
            end
            for i = 1:size(obj.oldtargets)
               obj.graphics(gcount) = plot(obj.oldtargets(i,1),obj.oldtargets(i,2),'r*');
               gcount = gcount + 1;
            end

            obj.graphics(gcount) = plot(obj.target(1),obj.target(2),'ys'); 
            % drawnow; % bottleneck  
        end
        
    end % methods
end % class

