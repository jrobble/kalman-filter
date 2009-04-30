% NOTE: When drawing rectangle, place person's face in center of upper half.

classdef hybrid_particle_set
    
    properties
        colorsd = 10; % standard deviation of color probability
        edgesd = 0.5; % standard deviation of edge probability
        origxvar = 10;
        origyvar = 10;

        orignumpts = 1000; % [100] number of particles in system
        
        % tracking variables
        graphics = []; % set of graphic handles
        steps = []; % number of steps to track
        sk = [];
        oldsk = [];
        oldtpos = [];
        oldnumpts = [];
        numpts = [];
        xvar = [];
        yvar = [];
        imgheight = [];
        imgwidth = [];
        islost = [];
        tpos = [];
        tcolormodel = [];
        tgxmodel = [];
        tgymodel = [];
        numedgebins = [];
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
        
        function obj = initialize(obj,img,iimg,numedgebins,igxchans,igychans)

            obj.imgheight = size(img,1);
            obj.imgwidth = size(img,2);
            obj.numedgebins = numedgebins;

            % select uniformly distributed random points w/in circle around target point
            fprintf(1,'Select ROI.\n\n');
            h = imrect(gca);
            fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
            setPositionConstraintFcn(h,fcn);
            obj.tpos = wait(h); % [x1,y1,width,height], [x1,y1] denotes upper left point
            delete(h);
            
            % calculate target model
            obj.tcolormodel = obj.calc_quad_feature_model(obj.tpos,iimg,3);
            obj.tgxmodel    = obj.calc_feature_model(obj.tpos,igxchans,obj.numedgebins);
            obj.tgymodel    = obj.calc_feature_model(obj.tpos,igychans,obj.numedgebins);
        end
        
        
        % step through alg.
        function [obj,islost] = track_target_step(obj,step,img,iimg,igxchans,igychans)
        
            % determine particle positions
            dist = obj.gen_uniform_norm_dist(obj.numpts,obj.tpos(1:2),obj.xvar,obj.yvar);
            
            % sk = [xcoord, ycoord, width, height, colorp, edgep]
            obj.sk = zeros(obj.numpts,6);
            obj.sk(:,1) = dist(:,1);
            obj.sk(:,2) = dist(:,2);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % STAGE 1 - color rectangles
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % initialize candidate pool
            candidates = 1:obj.numpts;
            
            % update candidates
            for i = candidates
                hpos(1:2) = obj.sk(i,1:2); % unique x and y
                hpos(3:4) = obj.tpos(3:4); % same width and height
                % hpos(3) = obj.tpos(3) + (obj.tpos(3) * (rand()/5 - 0.1));
                % hpos(4) = obj.tpos(4) + (obj.tpos(4) * (rand()/5 - 0.1));
                obj.sk(i,3:4) = hpos(3:4);
                
                % calculate color model
                hcolormodel = obj.calc_quad_feature_model(hpos,iimg,3);

                % handle color model
                colorp = -1;
                if ~isempty(hcolormodel) % particle not out of bounds
                   tk = obj.tcolormodel;
                   hk = hcolormodel;
                   rho = (tk-hk).^2;
                   rho = sum(rho,2);
                   rho = sum(rho);
                   rho = sqrt(rho);
                   colorp = exp(-(rho.^2)/obj.colorsd^2);
                end
                obj.sk(i,5) = colorp;
            end
            
            % determine probability distribution
            % keep top percentile and update candidates
            [colorpdist,colorpcenters] = hist(obj.sk(candidates,5),20); 
            [colorsortvals,colorsortindexes] = sort(obj.sk(candidates,5));
            colorsortindexes = candidates(colorsortindexes);
            
            [topindexes] = find(colorpdist); 
            topindex = topindexes(end);
            
            numcandidates = size(candidates,2);
            candidates = colorsortindexes((numcandidates-colorpdist(topindex)+1):numcandidates);
            % candidates % DEBUG
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % STAGE 2 - edge orientation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % update candidates
            for i = candidates
                hgxmodel = obj.calc_feature_model(hpos,igxchans,obj.numedgebins);
                hgymodel = obj.calc_feature_model(hpos,igychans,obj.numedgebins);

                % handle edge orientation model
                edgep = -1;
                if ~isempty(hgxmodel) && ~isempty(hgymodel) % particle not out of bounds
                    tk = [obj.tgxmodel,obj.tgymodel];
                    hk = [hgxmodel,hgymodel];
                    rho = (tk-hk).^2;
                    rho = sum(rho,2);
                    rho = sum(rho);
                    rho = sqrt(rho);
                    edgep = exp(-(rho.^2)/obj.edgesd^2);
                end
                obj.sk(i,6) = edgep;
            end
            
            % determine probability distribution
            % keep top percentile and update candidates
            [edgepdist,edgepcenters] = hist(obj.sk(candidates,6),10); 
            [edgesortvals,edgesortindexes] = sort(obj.sk(candidates,6));
            edgesortindexes = candidates(edgesortindexes);
            
            [topindexes] = find(edgepdist); 
            topindex = topindexes(end);
            
            numcandidates = size(candidates,2);
            candidates = edgesortindexes((numcandidates-edgepdist(topindex)+1):numcandidates);
            % candidates % DEBUG
            % pause();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % UPDATE TARGET
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % update target
            % [bestcolorp,bestcolorindex] = max(obj.sk(:,5));
            % [bestedgep, bestedgeindex]  = max(obj.sk(:,6));
            % bestindex = bestcolorindex;
            % bestp = bestcolorp;
            
            bestindex = candidates(end);
            bestp = obj.sk(bestindex,5);
            bestp
            
            %{
            fprintf(1,'bestcolorp(%d): %5.5f bestedgep(%d): %5.5f\n', ...
                bestcolorindex,bestcolorp,bestcolorindex,obj.sk(bestcolorindex,6));
            fprintf(1,'bestedgep(%d): %5.5f bestcolorp(%d): %5.5f\n', ...
                bestedgeindex,bestedgep,bestedgeindex,obj.sk(bestedgeindex,5));
            %}
            
            %{
            if bestp < 0.01 % assume target lost
                obj.islost = true;
            end
            %}
            
            obj.tpos(1:4) = obj.sk(bestindex,1:4);
            
            % newtk = obj.calc_quad_feature_model(obj.tpos,iimg,3);
            % obj.tk(1:2) = (bestp * newtk(1:2)) + ((1-bestp) * obj.tk(1:2));
            
            newtcolormodel = obj.calc_quad_feature_model(obj.tpos,iimg,3);
            newtgxmodel    = obj.calc_feature_model(obj.tpos,igxchans,obj.numedgebins);
            newtgymodel    = obj.calc_feature_model(obj.tpos,igychans,obj.numedgebins);
            
            obj.tcolormodel = (bestp * newtcolormodel) + ((1-bestp) * obj.tcolormodel);
            obj.tgxmodel    = (bestp * newtgxmodel)    + ((1-bestp) * obj.tgxmodel);
            obj.tgymodel    = (bestp * newtgymodel)    + ((1-bestp) * obj.tgymodel);
            
            
            % mark-up the image after determining points
            obj = obj.mark_up(step);
            islost = obj.islost;
        end
        
        
        function [fmodel] = calc_feature_model(obj,pos,iimg,numf)
            % one rectangular region
            width  = round(pos(3));
            height = round(pos(4));
            pos = round(pos);

            q1pos = [pos(1)+width,   pos(2)+height,   width, height];
            qpos = [q1pos];

            % check region validity
            for i = 1:size(qpos,1)
                if qpos(i,1) < 2 || qpos(i,2) < 2 || ...
                   qpos(i,1) - width > obj.imgwidth - 1 || qpos(i,2) - height > obj.imgheight - 1 
                    fmodel = [];
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

            q1sum = iimg(q1pos(2)-1,q1pos(1)-1,:);

            q1k = q1sum;

            % get rid of upper left data outside of ROI
            if pos(1) > 1 % not along left side
                q1k(1:numf) = q1k(1:numf) - iimg(q1pos(2)-1,        q1pos(1)-q1pos(3), 1:numf);
            end
            if pos(2) > 1 % not along upper side
                q1k(1:numf) = q1k(1:numf) - iimg(q1pos(2)-q1pos(4), q1pos(1)-1,        1:numf);
            end
            if pos(1) > 1 && pos(2) > 1
                q1k(1:numf) = q1k(1:numf) + iimg(q1pos(2)-q1pos(4), q1pos(1)-q1pos(3), 1:numf);
            end

            tmpq1k(1:numf) = q1k(:,:,1:numf)/(q1pos(3)*q1pos(4));

            q1k = tmpq1k;

            fmodel = [q1k];
        end
        
        
        function [fmodel] = calc_quad_feature_model(obj,pos,iimg,numf)
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
            for i = 1:size(qpos,1)
                if qpos(i,1) < 2 || qpos(i,2) < 2 || ...
                   qpos(i,1) - width > obj.imgwidth - 1 || qpos(i,2) - height > obj.imgheight - 1 
                    fmodel = [];
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
                q1k(1:numf) = q1k(1:numf) - iimg(q1pos(2)-1,        q1pos(1)-q1pos(3), 1:numf);
                q3k(1:numf) = q3k(1:numf) - iimg(q3pos(2)-1,        q3pos(1)-q3pos(3), 1:numf);
                q3k(1:numf) = q3k(1:numf) + iimg(q3pos(2)-q3pos(4), q3pos(1)-q3pos(3), 1:numf);
            end
            if pos(2) > 1 % not along upper side
                q1k(1:numf) = q1k(1:numf) - iimg(q1pos(2)-q1pos(4), q1pos(1)-1,        1:numf);
                q2k(1:numf) = q2k(1:numf) - iimg(q2pos(2)-q2pos(4), q2pos(1)-1,        1:numf);
                q2k(1:numf) = q2k(1:numf) + iimg(q2pos(2)-q2pos(4), q2pos(1)-q2pos(3), 1:numf);
            end
            if pos(1) > 1 && pos(2) > 1
                q1k(1:numf) = q1k(1:numf) + iimg(q1pos(2)-q1pos(4), q1pos(1)-q1pos(3), 1:numf);
            end

            tmpq1k(1:numf) = q1k(:,:,1:numf)/(q1pos(3)*q1pos(4));
            tmpq2k(1:numf) = q2k(:,:,1:numf)/(q2pos(3)*q2pos(4));
            tmpq3k(1:numf) = q3k(:,:,1:numf)/(q3pos(3)*q3pos(4));
            tmpq4k(1:numf) = q4k(:,:,1:numf)/(q4pos(3)*q4pos(4));

            q1k = tmpq1k;
            q2k = tmpq2k;
            q3k = tmpq3k;
            q4k = tmpq4k;

            fmodel = [q1k;q2k;q3k;q4k];
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

