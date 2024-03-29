% Taken from: http://www.mathworks.de/matlabcentral/fileexchange/17741
% By Ilya Rozenfeld
function[varargout]=draw_circle(rad, cntr, ltp, trnsp, img)

% DRAW_CIRCLE plots a circle patch which can be modified later by pressing 
% Shift+Left mouse button on the circle.
% 
% [CH, CNTR, RAD] = DRAW_CIRCLE(RAD, CNTR, LTP, TRNSP)
%
% Inputs:
%   RAD    -- radius
%   CNTR   -- 2-D aray of circle's center coordinates
%   LTP    -- color (default: blue)
%   TRNSP  -- degree of transparency  (default: transparent)
%
% All inputs are optional.  If input value needs to be skipped use [].  
% Radius and center coordinates are selected interactively by mouse 
% if not specified. Use left mouse button click to select postion and/or radius.
% In the interactive mode the information window will appear on the bottom
% of the figure window displaying relevant information about the circle
% patch.
%
% Outputs:
%   CH    --   handle to the created circle
%   CNTR  --   2-D aray of circle's center coordinates
%   RAD   --   radius
%
% When modification mode is entered (by pressing Shift+Left mouse button)
% modification dots (black squares) will appear.  The following operations 
% can be performed when in the modification mode:
%   -- The size and location of the circle can be changed by pressing on 
%      the dot and dragging the mouse.  Dots will disappear while the mouse
%      button is down and will reappear after mouse button is up.  
%   -- When pressing and holding left (right) mouse button the transparency 
%      of the circle patch will decrease (increase).  
%   -- Double-clicking on the patch will bring up color dialog box.
%
% While circle patch is being modified the information window will appear 
% on the bottom displaying relevant information.
% To exit modification mode click anywhere outside circle patch
%
% Created by Ilya Rozenfeld, 11/22/2007

CH = [];

imshow(img);
figProps = initialize_fig;

% Define defaults
if ~exist('ltp','var') || isempty(ltp)
    ltp = 'b';
end

if ~exist('trnsp','var') || isempty(trnsp)
    trnsp = 0;
end


% Select location of the circle's center.
if ~exist('cntr','var') || isempty(cntr)
    set(gcf, 'Pointer', 'fullcross')
    set(gcf, 'WindowButtonMotionFcn', @select_center)
    waitforbuttonpress;
    cp=get(gca,'CurrentPoint');
    x0=cp(1,1);
    y0=cp(1,2);
    set(gcf, 'Pointer', 'arrow')
else
    x0=cntr(1);
    y0=cntr(2);
end

% Select radius
if ~exist('rad','var') || isempty(rad)
    set(gcf, 'Pointer', 'fleur')
    set(gcf,'WindowButtonMotionFcn', @select_radius)
    waitforbuttonpress;
    set(gcf, 'Pointer', 'arrow')
else
    draw_fig(x0, y0, rad)
end

finalize_fig( figProps )

% Setup outputs
switch nargout

    case 1
        varargout{1} = CH;
    case 2
        varargout{1} = CH;
        varargout{2} = [x0,y0];
    case 3
        varargout{1} = CH;
        varargout{2} = [x0,y0];
        varargout{3} = rad;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function[] = modify_circle(src, eventdata)

        % Modify already plotted circle

        button = get(gcf, 'SelectionType');

        % The modification mode is activated by pressing
        % Shift+Left mouse button
        if ~strcmp(button,'extend') || ~isempty(findobj(gca,'Tag','modDots'))
            return
        end

        set(CH, 'ButtonDownFcn', '');
        figProps = initialize_fig;
        CH = gco;

        while 1

            % Plot modification dots
            h(1)=plot(x0,y0,'sk','markerfacecolor','k','markersize',8);
            h(2)=plot(x0,y0+rad,'sk','markerfacecolor','k','markersize',8);
            h(3)=plot(x0,y0-rad,'sk','markerfacecolor','k','markersize',8);
            h(4)=plot(x0+rad,y0,'sk','markerfacecolor','k','markersize',8);
            h(5)=plot(x0-rad,y0,'sk','markerfacecolor','k','markersize',8);

            x1=x0+rad*cos(pi/4);
            x2=x0-rad*cos(pi/4);
            y1=y0+rad*sin(pi/4);
            y2=y0-rad*sin(pi/4);
            h(6)=plot([x1 x1 x2 x2],[y2 y1 y1 y2],'sk','markerfacecolor','k','markersize',8);
            
            set(h,'Tag', 'modDots')

            try
                waitforbuttonpress;
            catch
                return
            end

            % Mouse click on the circle
            if gco == CH
                
                % Modify transparency and color of the object
                button = get(gcf, 'SelectionType');
                
                if strcmp(button, 'open')
                    
                    % Modify color
                    
                    circColor = uisetcolor(get(CH, 'EdgeColor'));
                    if length(circColor) > 1
                        set(CH, 'FaceColor', circColor)
                        set(CH, 'EdgeColor', circColor)
                    end
                    delete(h)
                    continue
                    
                elseif strcmp(button, 'normal')
                    transpFlg = 'inc';
                    
                elseif strcmp(button, 'alt')
                    transpFlg = 'dec';
                    
                else
                    delete(h)
                    continue
                    
                end

                setappdata(CH, 'TranspModify', 'on')
                button_up_fun = ...
                    @(src, eventdata) setappdata(CH, 'TranspModify', 'off');
                set(gcf, 'WindowButtonUpFcn', button_up_fun)
                delete(h)
                set_transp( transpFlg )
                set(gcf, 'WindowButtonUpFcn', '')
            
            elseif isempty(find(gco==h, 1))
                
                % Mouse click outside the cricle and modification dots
                % Exit modification mode
                finalize_fig( figProps )
                delete(h);
                break
                              
            else

                % Mouse click on one of the modification dots
                
                % Assign commands corresponding to the selected 
                % modification dot
                switch gco
                    case h(1)
                        set(gcf,'WindowButtonMotionFcn', @select_center)

                    case {h(2), h(3)}
                        set(gcf,'WindowButtonMotionFcn', @mod_radius_y)

                    case {h(4), h(5)}
                        set(gcf,'WindowButtonMotionFcn', @mod_radius_x)

                    case h(6)
                        set(gcf,'WindowButtonMotionFcn', @select_radius)
                end

                delete(h);

                set(gcf,'WindowButtonUpFcn', 'uiresume')
                uiwait
                set(gcf,'WindowButtonMotionFcn','')
                set(gcf,'WindowButtonUpFcn', '')
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function[ figSet ] = initialize_fig()

        % Setup figure porperties

        figSet.XMOD = get(gca,'xlimmode');
        figSet.YMOD = get(gca,'xlimmode');
        set(gca,'xlimmode','manual')
        set(gca,'ylimmode','manual')
        
        figSet.HLD = get(gca,'NextPlot');
        set(gca,'NextPlot','add')

        figSet.WBDF = get(gcf,'WindowButtonDownFcn');
        figSet.WBUF = get(gcf,'WindowButtonUpFcn');
        figSet.WBMF = get(gcf,'WindowButtonMotionFcn');
        set(gcf,'WindowButtonDownFcn','');
        set(gcf,'WindowButtonUpFcn','');
        set(gcf,'WindowButtonMotionFcn','');

        % Create message window
        figSet.HU=uicontrol(gcf,'Style','text',...
            'units','pixels',...
            'Position',[1 1 220 15],...
            'HorizontalAlignment','left');

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function[] = finalize_fig( figSet )

        % Restore figure properties after plotting is done

        delete(figSet.HU)
        set(gcf,'WindowButtonDownFcn', figSet.WBDF)
        set(gcf,'WindowButtonUpFcn', figSet.WBUF)
        set(gcf,'WindowButtonMotionFcn', figSet.WBMF)

        set(gca,'NextPlot', figSet.HLD)
        
        set(gca,'xlimmode', figSet.XMOD)
        set(gca,'ylimmode', figSet.YMOD)
        
        set(CH, 'ButtonDownFcn', @modify_circle);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function[] = select_center(src, eventdata)

        % Interactively select or modify circle's center
        cp=get(gca,'CurrentPoint');
        str=['X0 = ' num2str(cp(1,1),'%6.2f') ', Y0 = ' num2str(cp(1,2),'%6.2f')];
        set(figProps.HU,'String',str)

        if ~isempty(CH)
            x0 = cp(1,1);
            y0 = cp(1,2);
            draw_fig(x0, y0, rad)
        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function[] = select_radius(src, eventdata)

        % Interactively select or modify circle's radius

        cp=get(gca,'CurrentPoint');
        rad=sqrt((cp(1,1)-x0)^2+(cp(1,2)-y0)^2);
        str=['Radius = ' num2str(rad(1),'%6.2f')];
        set(figProps.HU,'String',str)
        draw_fig(x0, y0, rad)
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function[] = mod_radius_x(src, eventdata)

        % Interactively modify circle by dragging to the right or left

        cp=get(gca,'CurrentPoint');
        rad = (rad + abs(cp(1,1) - x0))/2;
        x0 = cp(1,1) - sign(cp(1,1) - x0)*rad;

        str=['X0 = ' num2str(x0,'%6.2f') ', Y0 = ' num2str(y0,'%6.2f')];
        str=[str ', Radius = ' num2str(rad,'%6.2f')];
        set(figProps.HU,'String',str)
        draw_fig(x0, y0, rad)
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function[] = mod_radius_y(src, eventdata)

        % Interactively modify circle by dragging up or down

        cp=get(gca,'CurrentPoint');
        rad = (rad + abs(cp(1,2) - y0))/2;
        y0 = cp(1,2) - sign(cp(1,2) - y0)*rad;

        str=['X0 = ' num2str(x0,'%6.2f') ', Y0 = ' num2str(y0,'%6.2f')];
        str=[str ', Radius = ' num2str(rad,'%6.2f')];
        set(figProps.HU,'String',str)
        draw_fig(x0, y0, rad)
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function[] = draw_fig(x0, y0, rad)
        
        % Draw circle patch

        t=0:0.01:2*pi;
        x=x0+rad*cos(t);
        y=y0+rad*sin(t);
        if isempty(CH)  % New patch
            CH = patch(x,y,ltp);
            set(CH,'facealpha',trnsp,'edgecolor',ltp);
            set(CH,'HitTest','off');
            
        else   % Modify patch
            set(CH,'xdata',x,'ydata',y);
        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function[] = set_transp(transpFlg)
        
        % Modify transparency of the patch
        
        currTransp = get(CH, 'FaceAlpha');
        transpInc = 0.01;
        maxTransp = 1;
        minTransp = 0;

        if strcmp(transpFlg,'inc')
            valRange = currTransp : transpInc : maxTransp;
        elseif strcmp(transpFlg, 'dec')
            valRange = currTransp: -transpInc : minTransp;
        end

        for transpVals = valRange
            if strcmp( getappdata(CH, 'TranspModify'), 'off')
                break
            end
            set(CH, 'FaceAlpha', transpVals)
            set(figProps.HU, 'String', ['Transparency = ' num2str(transpVals)])
            pause(0.1)
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end