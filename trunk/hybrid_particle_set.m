
classdef hybrid_particle_set
    
    properties
        sd = 10; % standard deviation of probability calc.
        origxvar = 10;
        origyvar = 10;

        orignumpts = 1000; % [100] number of particles in system
        
        % tracking variables
        graphics = []; % set of graphic handles
        steps = []; % number of steps to track
        sk = [];
        oldsk = [];
        oldtpos = [];
        target = [];
        oldnumpts = [];
        numpts = [];
        xvar = [];
        yvar = [];
        imgheight = [];
        imgwidth = [];
        islost = [];
        tpos = [];
        tk = [];
    end
    
    methods
        
        % constructor must be able to handle the no input case for object array pre-allocation
        function obj = hybrid_particle_set(steps)
            obj.islost = false;
            obj.numpts = obj.orignumpts;
            obj.oldnumpts = obj.numpts;
            obj.xvar = obj.origxvar;
            obj.yvar = obj.origyvar;
            if nargin > 0
                obj.steps = steps;
            end
        end
        
        function obj = initialize(obj,img,iimg)

            obj.imgheight = size(img,1);
            obj.imgwidth = size(img,2);

            % select uniformly distributed random points w/in circle around target point
            fprintf(1,'Select ROI.\n\n');
            h = imrect(gca);
            fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
            setPositionConstraintFcn(h,fcn);
            obj.tpos = wait(h); % [x1,y1,width,height], [x1,y1] denotes upper left point
            delete(h);
            
            % calculate target color model
            obj.tk = obj.calc_quad_color_model(obj.tpos,iimg);
        end
        
        
        % step through alg.
        function [obj,islost] = track_target_step(obj,step,img,iimg)
        
            % determine particle positions
            dist = obj.gen_uniform_norm_dist(obj.numpts,obj.tpos(1:2),obj.xvar,obj.yvar);
            
            % sk = [x,y,b,c] 
            % [xcoord, ycoord, probability, cumulative probability]
            obj.sk = zeros(obj.numpts,4);
            obj.sk(:,1) = dist(:,1);
            obj.sk(:,2) = dist(:,2);
            obj.sk(:,5) = 1;
            obj.sk(:,6) = 1;
            
            % update particles
            for i = 1:obj.numpts
                hpos(1:2) = obj.sk(i,1:2); % unique x and y
                hpos(3:4) = obj.tpos(3:4); % same width and height
                
                % sum color component diffs for each region
                hk = obj.calc_quad_color_model(hpos,iimg);
                
                if isempty(hk) % particle out of bounds
                   obj.sk(i,3) = -1; 
                else
                   rho = (obj.tk-hk).^2;
                   rho = sum(rho,2);
                   rho = sum(rho);
                   rho = sqrt(rho);
                   p = exp(-(rho.^2)/obj.sd^2);
                   obj.sk(i,3) = p;  
                end
                
                % fprintf(1,'pt: %d\n',i); obj.sk(i,3) % DEBUG
            end

            % update target
            [bestp,bestindex] = max(obj.sk(:,3));
            
            bestp % DEBUG
            if bestp < 0.01 % assume target lost
                obj.islost = true;
            end
            
            obj.tpos(1:2) = obj.sk(bestindex,1:2);
            obj.tk = obj.calc_quad_color_model(obj.tpos,iimg);
            
            %{
            tmpgraphics = [];
            for i = 1:obj.numpts
                % fprintf(1,'Calculating particle color distribution, [%d,%d] ...\n',obj.sk(i,1),obj.sk(i,2)); % DEBUG
                [py] = obj.get_color_distribution(obj.sk(i,:),img,obj.qccenters);

                if isnan(py) % DEBUG
                    obj.sk(i,:)
                    error('py is NAN');
                end
                
                % DEBUG
                
                fprintf(1,'obj.q: [ ');
                for j = 1:size(obj.q,2)
                    fprintf(1,'%1.2f ',obj.q(j));
                end
                fprintf(1,']\n');
                fprintf(1,'   py: [ ');
                for j = 1:size(py,2)
                    fprintf(1,'%1.2f ',py(j));
                end
                fprintf(1,']\n');
                
                %ro = 1 - sum(abs(py - obj.q))/2; % absolute error
                
                % calculate Bhattacharyya coefficient (ro)
                ro = sum(sqrt(py .* obj.q));
                fprintf(1,'rho: %d\n',ro);
                
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
                h = plot(obj.sk(i,1),obj.sk(i,2),'w+');
                g = obj.plot_circle(obj.rad,obj.sk(i,1),obj.sk(i,2));
                % obj.sk(i,3)
                pause();
                delete(h); 
                delete(g);
                % tmpgraphics(end+1) = h;
            end
            %}
            
            %pause();
            %for i = 1:size(tmpgraphics,2)
            %   delete(tmpgraphics(i)); 
            %end
            %}

            %{
            % normalize the probabilities
            unnormprobs = obj.sk(:,3);
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
            oldtarget = obj.target;
            if obj.targetmethod == 1
                % 1. "weighted mean" method 
                obj.target = [sum(obj.sk(:,3) .* obj.sk(:,1)), ...
                              sum(obj.sk(:,3) .* obj.sk(:,2))];
                obj.target % DEBUG
                
            elseif obj.targetmethod == 2
                % 2. "best particle" method
                [val,bestindex] = max(obj.sk(:,3));
                best(1:8) = obj.sk(bestindex,1:8);
                fprintf(1,'Best particle method best: %d\n',best(3)); % DEBUG
                obj.target(1:2) = best(1:2);
                
            elseif obj.targetmethod == 3
                % 3. "robust mean" method (weighted mean in small circular window around best particle)
                % plot(best(1),best(2),'w+'); % DEBUG
                [val,bestindex] = max(obj.sk(:,3));
                best(1:8) = obj.sk(bestindex,1:8);
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
            [obj.posqueue] = obj.posqueue.push_back(obj.target-oldtarget); % QUEUE
            
            % DEBUG - update new target color profile
            % if ~obj.islost
            %   oldq = obj.q;
            %    [obj.q,obj.qccenters] = obj.get_color_distribution(obj.target,img);
            %    obj.q = (oldq + obj.q) / 2;
            % end
                
            % DEBUG - mark points (particles) before update rule
            % xmin = min(obj.sk(:,1));
            % xmax = max(obj.sk(:,1));
            % ymin = min(obj.sk(:,2));
            % ymax = max(obj.sk(:,2));
            
            %{
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
            xdiff = obj.target(1) - oldtarget(1);
            ydiff = obj.target(2) - oldtarget(2);
            obj.scalex = xdiff;
            obj.scaley = ydiff;
            if obj.scalex < 5
                obj.scalex = 5 * sign(obj.scalex);
            end
            if obj.scaley < 5
                obj.scaley = 5 * sign(obj.scaley); 
            end

            maxb = best(3);
            if isnan(maxb) 
                error('maxb is NAN');
            end
            
            prob = unnormprobs(bestindex);
            fprintf(1,'unnormalized best: %d\n',prob); % DEBUG
            
            if prob < 0.70 % at least 70% similar
                if obj.islost == 0
                    % increase particle cloud size
                    fprintf(1,'>> Possible loss of target. ');
                    fprintf(1,'=> Use larger particle cloud and best particle method.\n');
                    obj.scalex = obj.origscalex * 1.5;
                    obj.scaley = obj.origscaley * 1.5;
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
                    obj.targetmethod = 3; % [3] use "robust mean" method
                    
                    % might want to scan the whole image here...
                    % reset Brownian motion parameters ?
                    % increase number of particles ?
                    obj.oldnumpts = obj.numpts;
                    % obj.numpts = 2 * obj.orignumpts;
                    obj.islost = 0;
                end
            end

            
            % select resample particles (can choose repeat particles)
            fprintf(1,'\nResample particles ...\n'); % DEBUG
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
            % sk1 = obj.sk; % DEBUG
            
            sk1 = zeros(obj.numpts,size(obj.sk,2));
            for i = 1:obj.numpts
                isvalid = 0;
                attempts = 0;
                while isvalid == 0
                    %{
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
                    sk1(i,1) = round((sk1(i,5)-sk0(i,5)) * obj.scalex + sk0(i,1));
                    sk1(i,2) = round((sk1(i,6)-sk0(i,6)) * obj.scaley + sk0(i,2));
                    %}
                    
                    % DEBUG - linear motion
                    % sk1(i,1) = sk0(i,1) + ((rand()-0.2) * obj.scalex * 2);
                    % sk1(i,2) = sk0(i,2) + ((rand()-0.2) * obj.scaley * 2); 
                    
                    [obj.posqueue,val] = obj.posqueue.mean(); % QUEUE
                    sk1(i,1) = sk0(i,1) + val(1) + (rand()-0.5) * 5;
                    sk1(i,2) = sk0(i,2) + val(2) + (rand()-0.5) * 5; 
                    
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
            % obj.oldsk = obj.sk;
            % obj.sk = sk1;
            %}
            
            % mark-up the image after determining points
            obj = obj.mark_up(step);
            islost = obj.islost;
        end
        
        
        function [cmodel] = calc_quad_color_model(obj,pos,iimg)
            % divide into equal-sized quadrants
            width  = round(pos(3)/2);
            height = round(pos(4)/2);
            pos = round(pos);
            
            q1pos = [pos(1)+width,   pos(2)+height,   width, height];
            q2pos = [pos(1)+2*width, pos(2)+height,   width, height];
            q3pos = [pos(1)+width,   pos(2)+2*height, width, height];
            q4pos = [pos(1)+2*width, pos(2)+2*height, width, height];
            qpos = [q1pos;q2pos;q3pos;q4pos];
            
            % check region validity
            for i = 1:size(qpos,2)
                if qpos(i,1) < 2 || qpos(i,2) < 2 || ...
                   qpos(i,1) - width > obj.imgwidth - 1 || qpos(i,2) - height > obj.imgheight - 1 
                    cmodel = [];
                    return;
                end
            end

            % update region bounds if necessary
            for i = 1:size(qpos,1)
                if qpos(i,1) > obj.imgwidth
                   qpos(i,2) = qpos(1,2) - (qpos(i,1) - obj.imgwidth); 
                   qpos(i,1) = obj.imgwidth;
                end
                if qpos(i,2) > obj.imgheight
                   qpos(i,4) = qpos(1,4) - (qpos(i,2) - obj.imgheight); 
                   qpos(i,2) = obj.imgheight; 
                end
            end
            
            q1pos(1:4) = qpos(1,1:4);
            q2pos(1:4) = qpos(2,1:4);
            q3pos(1:4) = qpos(3,1:4);
            q4pos(1:4) = qpos(4,1:4);
            
            q1sum = iimg(q1pos(2)-1,q1pos(1)-1,:);
            q2sum = iimg(q2pos(2)-1,q2pos(1)-1,:);
            q3sum = iimg(q3pos(2)-1,q3pos(1)-1,:);
            q4sum = iimg(q4pos(2)-1,q4pos(1)-1,:);
            
            q1k = q1sum;
            q2k = q2sum-q1sum;
            q3k = q3sum-q1sum;
            q4k = (q4sum+q1sum)-(q2sum+q3sum); % contains no excess
            
            % get rid of upper left data outside of ROI
            if pos(1) > 1 % not along left side
                q1k(1:3) = q1k(1:3) - iimg(q1pos(2)-1,        q1pos(1)-q1pos(3), 1:3);
                q3k(1:3) = q3k(1:3) - iimg(q3pos(2)-1,        q3pos(1)-q3pos(3), 1:3);
                q3k(1:3) = q3k(1:3) + iimg(q3pos(2)-q3pos(4), q3pos(1)-q3pos(3), 1:3);
            end
            if pos(2) > 1 % not along upper side
                q1k(1:3) = q1k(1:3) - iimg(q1pos(2)-q1pos(4), q1pos(1)-1,        1:3);
                q2k(1:3) = q2k(1:3) - iimg(q2pos(2)-q2pos(4), q2pos(1)-1,        1:3);
                q2k(1:3) = q2k(1:3) + iimg(q2pos(2)-q2pos(4), q2pos(1)-q2pos(3), 1:3);
            end
            if pos(1) > 1 && pos(2) > 1
                q1k(1:3) = q1k(1:3) + iimg(q1pos(2)-q1pos(4),q1pos(1)-q1pos(3),1:3);
            end
            
            tmpq1k(1:3) = q1k(:,:,1:3)/(q1pos(3)*q1pos(4));
            tmpq2k(1:3) = q2k(:,:,1:3)/(q2pos(3)*q2pos(4));
            tmpq3k(1:3) = q3k(:,:,1:3)/(q3pos(3)*q3pos(4));
            tmpq4k(1:3) = q4k(:,:,1:3)/(q4pos(3)*q4pos(4));
            
            q1k = tmpq1k;
            q2k = tmpq2k;
            q3k = tmpq3k;
            q4k = tmpq4k;
            
            cmodel = [q1k;q2k;q3k;q4k];
        end

        
        % TODO - may want to consider quasi-random sequence generator
        function [dist] = gen_uniform_norm_dist(obj,numpts,M,xvar,yvar)
            % M = [0 0];
            % C = [1 0; 0 1];
            C = [xvar 0; 0 yvar];
            dist = repmat(M,numpts,1) + randn(numpts,2)*C; % randn returns N by M matrix
            % figure; plot(dist(:,1),dist(:,2),'b+'); % DEBUG
        end
        

        %{
        % NOTE: function only used by first target
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
                
                % OLD WAY - use local ROI to determine centers
                %{
                area = floor(pi*obj.rad^2);
                innerpix = zeros(area,3);
                count = 1;
                for y = ymin:ymax
                    for x = xmin:xmax
                        % determine if pixel within circle ROI
                        if (x-centerx)^2+(y-centery)^2 <= obj.rad^2
                            innerpix(count,:) = img(y,x,:);
                            count = count + 1;
                        end

                    end
                end
                count = count - 1;
                innerpix = innerpix(1:count,:);
                % fprintf(1,'innerpix: %d\n',count); % DEBUG

                % calculate color histogram centers for ROI
                [rdist,rcenters] = hist(innerpix(:,1),obj.numbins);
                [gdist,gcenters] = hist(innerpix(:,2),obj.numbins);
                [bdist,bcenters] = hist(innerpix(:,3),obj.numbins);
                ccenters = [rcenters;gcenters;bcenters];
                %}
                
                [rdist,rcenters] = hist(0:255,obj.numbins);
                [gdist,gcenters] = hist(0:255,obj.numbins);
                [bdist,bcenters] = hist(0:255,obj.numbins);
                ccenters = [rcenters;gcenters;bcenters];
                ccenters % DEBUG
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

        
        % NOTE: function only used by first target
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
            
            if obj.usehsv == 1
                cprofile = [posr-1,posg-1,0];
            else
                cprofile = [posr-1,posg-1,posb-1];
            end
            
            % cprofile % DEBUG
        end
        %}
        
        % mark-up the image
        function obj = mark_up(obj,step)
           % clear previous graphics
            for i = 1:size(obj.graphics,2)
               delete(obj.graphics(i)); 
            end
            
            % mark points (particles)
            gcount = 1;
            % obj.graphics(gcount) = plot(obj.oldsk(:,1),obj.oldsk(:,2),'b+');
            % gcount = gcount + 1;
            
            % obj.graphics(gcount) = plot(obj.sk(:,1),obj.sk(:,2),'g+');
            % gcount = gcount + 1;

            % mark targets (tracking)
            % if ~obj.islost % only record good targets
                obj.oldtpos(end+1,:) = obj.tpos(:);
            % end
            for i = 1:size(obj.oldtpos)
               obj.graphics(gcount) = rectangle('Position',obj.tpos,'LineWidth',1,'EdgeColor','r'); % DEBUG
               gcount = gcount + 1;
            end

            obj.graphics(gcount) = rectangle('Position',obj.tpos,'LineWidth',1,'EdgeColor','y'); % DEBUG 
            % drawnow; % bottleneck  
        end
        
    end % methods
end % class

