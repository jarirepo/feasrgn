classdef feasrgnplot < handle
% feasrgnplot Feasible region plot
% 
% feasrgnplot Properties:
%   fillcolor
%   facealpha
%   linestyle
%   linewidth
%   edgecolor
%   nodelabel
%   nodeshape
%   nodesize
%   nodecolor
%   nodevalueformat
%   displaynodes
%   displaynodevalues
%   displayconstraints
%   displaylegend
%   displayinnerpoints
%   upointcount
%   vpointcount
% 
% Features:
%   -Allows direct manipulation of the plot style properties
%{
    Version : 1.0
    Date    : 2016-12-08
    Author  : Jari Repo, University West, jari.repo@hv.se
%}
    properties (SetObservable, AbortSet = true)
        fillcolor = [0 .75 .15]
        facealpha = 1
        linestyle = '-'
        linewidth = 1
        edgecolor = [0 0 0]
        nodelabel = 'S'
        nodeshape = 'o'
        nodesize = 5
        nodecolor = [0 0 .75]
        nodefontsize = 10
        nodevalueformat = '(%.2f; %.2f)'
        displaynodes = true
        displaynodevalues = true
        displayconstraints = true
        displaylegend = true;       
        displayquerypoint = true;
        displayinnerpoints = false;
        upointcount = 30;
        vpointcount = 30;
    end    
    properties (SetAccess = private, Hidden = true)
        handles = struct.empty
        frgn = feasrgn.feasrgn.empty;
    end
    properties (Constant)
        MAX_PTS = 50
    end
    methods
        function obj = feasrgnplot(frgn, varargin)
            % frgn - feasrgn object
            if nargin>0
                if isa(frgn,'feasrgn.feasrgn')
                    if frgn.closed
                        obj.frgn = frgn;
                        hold on
                        obj.handles = struct();
                        % output constraints functions
                        obj.handles.ph_cfcn = plot(frgn.x,frgn.F,'-',...
                            'LineWidth',obj.linewidth,...
                            'Visible','on');
                        nfcn = length(frgn.ineq);
                        for i=1:nfcn
                            set(obj.handles.ph_cfcn(i),'LineStyle',...
                                frgn.ineq(i).getLineStyle())
                        end
                        % output figure legend
                        legStr = cell(1,nfcn);
                        for i=1:nfcn
                            % get function signature for figure legend                
                            legStr{i} = frgn.ineq(i).getSignature();
                        end                                    
                        obj.handles.leg = legend(legStr,...
                            'FontSize',get(gca,'FontSize'),...
                            'Box','on',...
                            'Location','eastoutside',...
                            'Visible','on');
%                         set(get(obj.handles.leg,'String'),'FontSize',get(gca,'FontSize'))
                        % output a filled feasible region
                        % (area sometimes produces weird results)
%                         obj.handles.hfrgn = area(frgn.B(1,1:end-1),frgn.B(2,1:end-1),...
%                             'FaceColor',obj.fillcolor,...
%                             'LineStyle','none',...
%                             'LineWidth',obj.linewidth,...
%                             'EdgeColor','none',...
%                             'Visible','on');
                        Nb = size(frgn.B,2);
                        obj.handles.hfrgn = patch(frgn.B(1,:),frgn.B(2,:),-1e-6*ones(1,Nb),...
                            obj.fillcolor,...
                            'Clipping','on',...
                            'LineStyle',obj.linestyle,...
                            'LineWidth',obj.linewidth,...
                            'EdgeColor',obj.edgecolor,...
                            'FaceAlpha',obj.facealpha,...
                            'Visible','on');                        
                        % output feasible region boundary
                        if ~isempty(frgn.B)
                            obj.handles.ph_B = plot(frgn.B(1,:),frgn.B(2,:),...
                                obj.linestyle,...
                                'LineWidth',obj.linewidth,...
                                'Color',obj.edgecolor,...
                                'Visible','off');
                        end
                        if ~isempty(frgn.S)
                            % output nodes
                            obj.handles.ph_S = plot(frgn.S(1,:),frgn.S(2,:),...
                                obj.nodeshape,...
                                'MarkerSize',obj.nodesize,...
                                'Color',obj.nodecolor,...
                                'MarkerFaceColor',obj.nodecolor,...
                                'Visible','on');
                            % output node labels and values
                            ns = size(frgn.S,2);
                            str = cell(1,ns);
                            obj.handles.th_S = zeros(1,ns);
                            for i=1:ns
                                s = sprintf('%s_{%d}=%s',...
                                    obj.nodelabel, i, ...
                                    obj.nodevalueformat);
                                str{i} = sprintf(s,frgn.S(1,i), frgn.S(2,i));
                                obj.handles.th_S(i) = text(frgn.S(1,i),frgn.S(2,i),...
                                    str{i},...
                                    'FontSize',obj.nodefontsize,...
                                    'FontWeight','normal',...
                                    'HorizontalAlignment','left',...
                                    'VerticalAlignment','bottom',...
                                    'Clipping','on',...
                                    'Visible','on');
                            end                            
                        end                        
                        % output query point (initially hidden)
                        obj.handles.ph_Q = plot3(0,0,1e-6,'d',...
                            'MarkerSize',6,...
                            'Color',[1 .5 .5],...
                            'MarkerFaceColor',[1 .5 .5],...
                            'Visible','off');
                        % output inner points (initially hidden)
                        obj.handles.ph_Pts = plot3(0,0,0,'o',...
                            'MarkerSize',1,...
                            'Color',[0 0 .25],...
                            'MarkerFaceColor',[0 0 .25],...
                            'Visible','off');
                        % set appropriate Zoom into the feasible region
                        xmin = min(frgn.B(1,:)); xmax = max(frgn.B(1,:));
                        ymin = min(frgn.B(2,:)); ymax = max(frgn.B(2,:));
                        dx = (xmax-xmin)/10;
                        dy = (ymax-ymin)/10;
                        axis([xmin-dx xmax+dx ymin-dy ymax+dy])
                        hold off
                    else
                        error('A closed feasible region boundary is required')
                    end
                else
                    error('Input FRGN must be a feasrgn object')
                end
            end
            % add property listener            
            addlistener(obj,...
                {'fillcolor','facealpha','linestyle','linewidth','edgecolor','nodelabel',...
                'nodeshape','nodesize','nodecolor','nodefontsize','displaynodes','nodevalueformat',...
                'displaynodevalues','displayconstraints','displaylegend','displayquerypoint',...
                'displayinnerpoints','upointcount','vpointcount'},...
                'PostSet',@(src,evt)feasrgn.feasrgnplot.propertyChangeListener(obj,src,evt,obj.handles));
            % call property set methods
            for k=1:2:length(varargin)-1
                obj.(lower(varargin{k})) = varargin{k+1};
            end
            % add query point listener
            addlistener(frgn,'queryPt','PostSet',@(src,evt)feasrgn.feasrgnplot.propertyChangeListener(obj,src,evt,obj.handles));
        end
        
        function set.fillcolor(obj, value)
            obj.fillcolor = value;
        end
        function set.facealpha(obj, value)
            if isnumeric(value) && isreal(value) && value>=0 && value<=1            
                obj.facealpha = value;
            else
                error('Invalid value')
            end
        end
        function set.linestyle(obj, value)
            obj.linestyle = value;
        end
        function set.linewidth(obj, value)
            obj.linewidth = value;
        end
        function set.edgecolor(obj, value)
            obj.edgecolor = value;
        end
        function set.nodelabel(obj, value)
            if isa(value,'char') && ~isempty(strtrim(value))
                obj.nodelabel = value;
            else
                error('Invalid value')
            end
        end
        function set.nodeshape(obj, value)
            switch value
                case {'.','o','x','+','*','s','d','v','^','<','>','p','h'}
                    obj.nodeshape = value;
                otherwise
                    error('Invalid value')                    
            end
        end
        function set.nodesize(obj, value)
            if isnumeric(value) && isreal(value) && value>0
                obj.nodesize = value;
            else
                error('Invalid value')                    
            end
        end
        function set.nodecolor(obj, value)            
            obj.nodecolor = value;
        end
        function set.nodefontsize(obj, value)            
            obj.nodefontsize = value;
        end
        function set.displaynodes(obj, value)
            if islogical(value)
                obj.displaynodes = value;
            else
                error('Invalid value')                    
            end
        end
        function set.displaynodevalues(obj, value)
            if islogical(value)
                obj.displaynodevalues = value;
            else
                error('Invalid value')                    
            end
        end
        function set.displayconstraints(obj, value)
            if islogical(value)
                obj.displayconstraints = value;
            else
                error('Invalid value')                    
            end
        end
        function set.displaylegend(obj, value)
            if islogical(value)
                obj.displaylegend = value;
            else
                error('Invalid value')                    
            end
        end
        function set.displayquerypoint(obj, value)
            if islogical(value)
                obj.displayquerypoint = value;
            else
                error('Invalid value')                    
            end
        end
        function set.nodevalueformat(obj, value)
            if ischar(value) && ~isempty(value)
                obj.nodevalueformat = value;
            else
                error('Invalid value')                    
            end
        end        
        function set.displayinnerpoints(obj, value)
            if islogical(value)
                obj.displayinnerpoints = value;
            else
                error('Invalid value')                    
            end
        end
        function set.upointcount(obj, value)
            if feasrgn.feasrgnplot.islimint(value,2,feasrgn.feasrgnplot.MAX_PTS)
                obj.upointcount = value;
            else
                error('Invalid value')                    
            end
        end
        function set.vpointcount(obj, value)
            if feasrgn.feasrgnplot.islimint(value,2,feasrgn.feasrgnplot.MAX_PTS)
                obj.vpointcount = value;
            else
                error('Invalid value')                    
            end
        end
    end
        
    methods (Static)
        function value = islimint(a, amin, amax)
            value = ( isnumeric(a) && isreal(a) && numel(a)==1 && a==fix(a) && a>=amin && a<=amax );
        end
        function propertyChangeListener(hObj, hSrc, evtData, handles)
            switch hSrc.Name
                case 'queryPt'                    
                    set(handles.ph_Q,...
                        'XData',hObj.frgn.queryPt(1),...
                        'YData',hObj.frgn.queryPt(2))                                    
                case 'fillcolor'
                    set(handles.hfrgn,'FaceColor',hObj.fillcolor)                    
                case 'facealpha'
                    set(handles.hfrgn,'FaceAlpha',hObj.facealpha)                    
                case 'linestyle'
                    set(handles.hfrgn,'LineStyle',hObj.linestyle)
                case 'linewidth'
                    set(handles.hfrgn,'LineWidth',hObj.linewidth)
                case 'edgecolor'
                    set(handles.hfrgn,'EdgeColor',hObj.edgecolor)                    
                case 'nodelabel'
                    % update the node labels
                    for h=handles.th_S
                        str = get(h,'String');
                        res = regexp(str,'([a-zA-Z]*)[^_]','match','once'); 
                        strNew = strrep(str,res,hObj.nodelabel);
                        set(h,'String',strNew)
                    end
                case 'nodevalueformat'
                    % update the node values
                    for i=1:length(handles.th_S)
                        s = sprintf('%s_{%d}=%s',...
                            hObj.nodelabel, i, ...
                            hObj.nodevalueformat);
                        str = sprintf(s, hObj.frgn.S(1,i), hObj.frgn.S(2,i));                        
                        set(handles.th_S(i),'String',str)
                    end                    
                case 'nodeshape'
                    set(handles.ph_S,'Marker',hObj.nodeshape)
                case 'nodesize'
                    set(handles.ph_S,'MarkerSize',hObj.nodesize)
                case 'nodecolor'
                    set(handles.ph_S,'Color',hObj.nodecolor,...
                        'MarkerFaceColor',hObj.nodecolor)
                case 'nodefontsize'
                    set(handles.th_S,'FontSize',hObj.nodefontsize)
                case 'displaynodes'
                    if hObj.displaynodes
                        set(handles.ph_S,'Visible','on')
%                         % show node values
%                         if hObj.displaynodevalues
%                             set(handles.th_S,'Visible','on')
%                         end
                    else
                        set(handles.ph_S,'Visible','off')
%                         % hide node values as well
%                         set(handles.th_S,'Visible','off')
                    end
                case 'displaynodevalues'
                    if hObj.displaynodevalues
                        set(handles.th_S,'Visible','on')
                    else
                        set(handles.th_S,'Visible','off')
                    end
                case 'displayconstraints'
                    if hObj.displayconstraints
                        set(handles.ph_cfcn,'Visible','on')
                    else
                        set(handles.ph_cfcn,'Visible','off')
                    end
                case 'displaylegend'
                    if hObj.displaylegend
                        set(handles.leg,'Visible','on')
                    else
                        set(handles.leg,'Visible','off')
                    end
                case 'displayquerypoint'
                    if hObj.displayquerypoint
                        set(handles.ph_Q,'Visible','on')
                    else
                        set(handles.ph_Q,'Visible','off')
                    end
                case 'displayinnerpoints'
                    if hObj.displayinnerpoints
                        if isempty(hObj.frgn.Gx)
                            hObj.frgn.generateInnerPts(hObj.upointcount-1, hObj.vpointcount-1);
                            set(handles.ph_Pts,...
                                'XData',hObj.frgn.Gx(:),...
                                'YData',hObj.frgn.Gy(:),...
                                'ZData',2e-6*ones(size(hObj.frgn.Gx(:))))
                        end
                        set(handles.ph_Pts,'Visible','on')
                    else
                        set(handles.ph_Pts,'Visible','off')
                    end
                case {'upointcount','vpointcount'}
                    hObj.frgn.generateInnerPts(hObj.upointcount-1, hObj.vpointcount-1);
                    set(handles.ph_Pts,...
                        'XData',hObj.frgn.Gx(:),...
                        'YData',hObj.frgn.Gy(:),...
                        'ZData',2e-6*ones(size(hObj.frgn.Gx(:))))
            end
        end
    end
end
