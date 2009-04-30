% NOTE: When drawing rectangle, place person's face in center of upper half.

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
                % obj.islost = true;
            end
            
            obj.tpos(1:2) = obj.sk(bestindex,1:2);
            newtk = obj.calc_quad_color_model(obj.tpos,iimg);
            obj.tk = (bestp * newtk) + ((1-bestp) * obj.tk);

            
            
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
            
            obj.graphics(gcount) = plot(obj.sk(:,1),obj.sk(:,2),'g+');
            gcount = gcount + 1;

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

