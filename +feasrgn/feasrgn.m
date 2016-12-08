
classdef feasrgn < handle
% feasrgn Piecewise linear feasible region based on a set of inequalities
% 
% USAGE: frgn = feasrgn(cfcn, opts); 
%
%       cfcn - a cell array defining the inequalities
%       opts - an options struct with fields 'xmin','xmax','dx'
% 
% Example:
%     cfcn = { 'y','>', @(x) .65*(x-1).^2-3; ...
%              'y','<=',@(x) -1.25*x.^2+6; ...
%              'y','>=', @(x) .8*(x+1.25).^3-1.025*(x+.5); ...
%              'y','>', @(x) .5*ones(size(x)); ...
%              'y','<=',@(x) 3.5*ones(size(x)); ...
%              'y','<=',@(x) -3*x+2; ...
%              'y','<=',@(x) 3-.5*x };
%     opts = struct('xmin',-3,'xmax',3,'dx',0.05);
%     frgn = feasrgn(cfcn,opts);
% 
% feasrgn Properties (read-only)
%   ineq
%   x
%   F
%   S
%   B
%   closed
%   queryPt
% 
% feasrgn Methods:
%   isPtInside
%   generateInnerPts
% 
% Dependencies:
%   feasrgn.inequality
%   feasrgn.ineqsign

%{
    --- Change log ---

    Version 1.2, JRE
    -The inequalities now need to be defined in a cell array which allows
    feasrgn to automatically number the inequalities according to their 
    order of occurence in the cell array
    -Improved the scanning algorithm further by utilizing node connectivity
    information to ensure that all intersections will be included in the 
    final generated boundary
    -Colinear boundary segments are automatically removed
    -The cell-array with legend string (legStr) is created in feasrgnplot
    to keep this class clean from the visualization
    -Improved validation of the input
    -New method isPtInside to test if a query point is inside the feasible
    region boundary

    Version 1.1, JRE
    -Rewrote major portions of the scanning loop which resolved some issues
    with missing intersections in the generated boundary

    Version : 1.2
    Date    : 2016-12-06
    Author  : Jari Repo, University West, jari.repo@hv.se
%}

    properties (SetAccess = private) 
        ineq = feasrgn.inequality.empty % array of inequality objects
        closed = false      % indicates if B represents a closed boundary
        x = []              % discrete grid for the independent variable
        F = []              % evaluated constraints functions
        S = []              % intersections/nodes [x1 x2 ... xn; y1 y2 ... yn]
        B = []              % feasible region boundary [x1 x2 ... xm; y1 y2 ... ym]
        Gx = [];            % generated inner points (Gx,Gy)
        Gy = [];            %               "
    end    
    properties (SetObservable, AbortSet = true)
        queryPt = [0;0]
    end    
    properties (SetAccess = private, Hidden = true)
        p_ = []
        q_ = []
    end
    properties (Constant, Hidden = true)
        DEBUG = 0
        MAX_NODES = 99 
        COLIN_TOL = 1e-11
    end            
    methods
        function obj = feasrgn(cfcn, opts)
            % cfcn - a cell array defining the inequalities
            % opts - an options struct with fields 'xmin','xmax','dx'
            
            % check on the inputs
            if nargin>0
                if ~iscell(cfcn)
                    error('Input CFCN must be a cell array with constraints functions definitions')
                end
                [nfcn,nc] = size(cfcn);                
                if nc~=3
                    error('Input CFCN must be a 3-column cell array')
                end
                if nfcn<2
                    error('At least two constraints functions are required')
                end
                
                % validate cell array cfcn
                for i=1:nfcn
                    % 1st column
                    if ~isa(cfcn{i,1},'char')
                        error('Constraints function #%d: dVarName must be a char',i)
                    else
                        if isempty(strtrim(cfcn{i,1}))
                            error('Constraints function #%d: dVarName must be a char',i)
                        end
                    end
                    % 2nd column
                    if ~isa(cfcn{i,2},'char')
                        error('Constraints function #%d: inequality sign must be a char',i)                        
                    else
                        if isempty(strtrim(cfcn{i,2}))
                            error('Constraints function #%d: inequality sign must be a char',i)                        
                        else
                            % test enumeration
                            o = feasrgn.ineqsign.enum(cfcn{i,2});                            
                        end
                    end
                    % 3rd column
                    if ~isa(cfcn{i,3},'function_handle')
                        error('Constraints function #%d: hFcn must be a function_handle',i)
                    end
                end
                
                % validate options
                if nargin>1 && isstruct(opts)
                    if isfield(opts,'xmin')
                        xmin = opts.xmin;
                    else                        
                        error('Inputs OPTS must contain the field ''xmin''')
                    end
                    if isfield(opts,'xmax')
                        xmax = opts.xmax;
                    else                        
                        error('Inputs OPTS must contain the field ''xmax''')
                    end
                    if isfield(opts,'dx')
                        dx = opts.dx;
                    else                        
                        error('Inputs OPTS must contain the field ''dx''')
                    end
                    if xmin>=xmax
                        error('Invalid interval [xmin,xmax]')
                    end                    
                    r = xmax-xmin;
                    if dx<=0 || dx>=r
                        error('Invalid resolution dx specified')
                    end
                    x = (xmin : dx : xmax)';
                    
                else
                    % no options specified. Find out appropriate values for 
                    % [xmin,xmax] by using a "hidden ezplot"
                    xmin = Inf;
                    xmax = -Inf;
                    dx = Inf;                    
                    fh = figure('Visible','off');
                    for i=1:nfcn
                        if ~isa(cfcn{i,3},'function_handle')
                            error('Missing function handle for constraints function #%d',i)
                        end
                        x = get( ezplot(cfcn{i,3}),'XData');
                        x1 = min(x);
                        x2 = max(x);
                        if x1<xmin, xmin=x1; end
                        if x2>xmax, xmax=x2; end
                        h = x(2)-x(1);
                        if h<dx, dx=h; end
                    end
                    close(fh);
                    % increase resolution 10 times
                    dx = dx/10;
                    x = (xmin : dx : xmax)';
                    xmax = x(end);
                    warning(sprintf('Using xmin=%.4f, xmax=%.4f, dx=%.4f\n',...
                        xmin,xmax,dx))
                end
            else
                error('Missing inputs')
            end
            
            % create numbered inequality objects            
            %   -the numbering is based on their order of occurence in the 
            %   cfcn cell array
            N = length(x);
            F = zeros(N,nfcn);
            for i=1:nfcn
                obj.ineq(i) = feasrgn.inequality(i,cfcn{i,1},cfcn{i,2},cfcn{i,3});
                % evaluate constraints function
                F(:,i) = obj.ineq(i).hFcn(x);
            end
            
            % sort the constraints functions in ascending order
            [Fs,Isrt] = sort(F,2,'ascend');
            
            % scan for intersections...
            S = [];
            ns = 0;
            Fx = zeros(1,nfcn);
            dtol = 1e-11;
            X = x;
            Fupd = F;
            nodes(1:feasrgn.feasrgn.MAX_NODES) = struct('conn',[]);                        
%             % init level nodes
%             levels(1:nfcn) = struct('level',NaN,'nodes',[]);
%             for i=1:nfcn
%                 levels(i).level = i;                
%             end
            
            for k=1:N-1
                I = Isrt(k,:);
                for i=1:nfcn-1
                    a = I(i);
                    Pa = [x(k),x(k+1); F(k,a),F(k+1,a)];                
                    for j=i+1:nfcn
                        b = I(j);
%                         if a==b,error('Unexpected error'),end
                        Pb = [x(k),x(k+1); F(k,b),F(k+1,b)];
                        % intersect functions a and b
                        A = [Pa(:,2)-Pa(:,1), -(Pb(:,2)-Pb(:,1))];            
                        if abs(det(A))>1e-11
                            t = A\(Pb(:,1)-Pa(:,1));
                            if t(1)>-1e-11 && t(1)<(1+1e-11)
                                % intersection found
                                Q = Pa(:,1)+(Pa(:,2)-Pa(:,1))*t(1);
                                % avoid dublicate intersections                                
                                nodeIndex = [];
                                if ~isempty(S)
                                    d = sqrt(sum((Q*ones(1,ns)-S).^2));
                                    nodeIndex = find(d < dtol);
%                                     if length(nodeIndex)>1,error('Unexpected error'),end
                                end       
                                if isempty(nodeIndex)
                                    % evaluate the CF:s at the intersection
                                    for m=1:nfcn
                                        if m==a || m==b
                                            Fx(m) = Q(2);
                                        else
                                            Fx(m) = obj.ineq(m).hFcn(Q(1));
                                        end
                                    end   
                                    % save unique node Q
                                    S = [S, Q];
                                    ns = ns+1;                                    
                                    X = [X; Q(1)];
                                    Fupd = [Fupd; Fx];  
                                    % register inequalities a and b to the new node
                                    if ~any(nodes(ns).conn==a)
                                        nodes(ns).conn = [nodes(ns).conn, a];
                                    end
                                    if ~any(nodes(ns).conn==b)
                                        nodes(ns).conn = [nodes(ns).conn, b];
                                    end                                                                        
                                else
                                    % node already exists
                                    % register inequalities a and b to the
                                    % existing node
                                    if ~any(nodes(nodeIndex).conn==a)
                                        nodes(nodeIndex).conn = [nodes(nodeIndex).conn, a];
                                    end
                                    if ~any(nodes(nodeIndex).conn==b)
                                        nodes(nodeIndex).conn = [nodes(nodeIndex).conn, b];
                                    end                                    
                                end
                            end
                        end
                    end
                end
            end
            
            if feasrgn.feasrgn.DEBUG
                fprintf(1,'Number of nodes: %d\n',ns);
                for i=1:ns
                    fprintf(1,'Node S%d: ',i);
                    for j=nodes(i).conn
                        fprintf(1,'%6d',j);
                    end
                    fprintf(1,'\n');
                end
            end            
            
            
            [Xs,I] = sort(X,'ascend');
            F = Fupd(I,:);
            [Fs,Isrt] = sort(F,2,'ascend');
            
            % get indices for the intersections
%             isectIdx = find(I>N)
            
            obj.x = Xs;
            obj.F = F;
            obj.S = S;            
            
            Bx = [];
            By = [];
            Bsgn = [];
            Bfcn = [];
            
            % scan the sorted data array for the feasible region...
            if feasrgn.feasrgn.DEBUG
                fprintf(1,'Scanning in progress...\n');
            end
            sgn = zeros(1,nfcn);
            N = length(Xs);
            validx = zeros(1,N); validx(:) = false;
            iLow = zeros(1,N);
            
            for k=1:N
                % ignore unbounded regions "downwards" or "upwards"    
                % (potential intersection will be missed but fixed in the
                % next iteration)
                if validx(k) || obj.ineq(Isrt(k,1)).isNeg() || obj.ineq(Isrt(k,nfcn)).isPos()
                    continue
                end
                % convert ineq. sign into a numerical value (-1 or +1)                
                nChanges = 0;
                for i=1:nfcn
                    sgn(i) = obj.ineq(Isrt(k,i)).ineqSign.getDirection();
                    if i>1
                        if sgn(i-1)*sgn(i)<0
%                         if sgn(i-1)==1 && sgn(i)==-1
                            nChanges = nChanges+1;
                            if nChanges==1
                                % record the index for the first occurence of a sign change
                                % only one sign change is allowed for a feasible region                            
                                iLow(k) = i-1;
                            else
                                % exit this for-loop as soon as more than one
                                % sign change is encountered
                                break;
                            end
                        end
                    end
                end                
                if nChanges==1           
                    % check if the previous X was a missed intersection/node
                    % using the node connectivity information
                    if k>1 && ns>0 && ~validx(k-1)
                        % get previous point on this CF
                        P = [Xs(k-1); Fs(k-1,iLow(k))];
                        nodeFound = false;
                        for i=1:ns
                            for j=nodes(i).conn
                                if j==Isrt(k,iLow(k))
                                    d = sqrt(sum((P-S(:,i)).^2));
                                    if d<dtol                                        
                                        Bx = [Bx,S(1,i)];
                                        By = [By,S(2,i)];
                                        Bsgn = [Bsgn,sgn(iLow(k))];
                                        Bfcn = [Bfcn,Isrt(k,iLow(k))];
                                        validx(k-1) = true;
                                        iLow(k-1) = iLow(k);                                        
                                        nodeFound = true;
                                        break
                                    end
                                end
                            end
                            if nodeFound,break,end
                        end
                    end
                    
                    % add new boundary point (avoid multiples)
                    d = Inf;
                    if numel(Bx)>0
                        d = sqrt((Xs(k)-Bx(end))^2+(Fs(k,iLow(k))-By(end))^2);
                    end
                    if d>dtol
                        Bx = [Bx,Xs(k)];
                        By = [By,Fs(k,iLow(k))];
                        Bsgn = [Bsgn,sgn(iLow(k))];
                        Bfcn = [Bfcn,Isrt(k,iLow(k))];                        
                        validx(k) = true;
                    end
                    
                    % look one step ahead so we don't miss any later
                    % intersection when traversing this CF    
                    if k<N && ns>0
                        % get next point on this CF
                        P = [Xs(k+1); Fs(k+1,iLow(k))];
                        nodeFound = false;
                        for i=1:ns
                            for j=nodes(i).conn
                                if j==Isrt(k,iLow(k))
                                    d = sqrt(sum((P-S(:,i)).^2));
                                    if d<dtol                                        
                                        Bx = [Bx,S(1,i)];
                                        By = [By,S(2,i)];
                                        Bsgn = [Bsgn,sgn(iLow(k))];
                                        Bfcn = [Bfcn,Isrt(k,iLow(k))];                                        
                                        validx(k+1) = true;
                                        iLow(k+1) = iLow(k);                                        
                                        nodeFound = true;
                                        break
                                    end
                                end
                            end
                            if nodeFound,break,end
                        end                        
                    end
                end                
            end                        

            % generate the lower and upper boundaries for the ruled surface
            J = find(validx);            
            if ~isempty(J)
                xx = Xs(J);
                p = Fs(J,iLow(J));
                obj.p_ = feasrgn.feasrgn.removeColinear([xx,p(:,1)]');
                q = Fs(J,iLow(J)+1);
                obj.q_ = feasrgn.feasrgn.removeColinear([xx,q(:,1)]');
                %figure,plot(obj.p_(1,:),obj.p_(2,:),'.-',obj.q_(1,:),obj.q_(2,:),'.-')
            end
            
            % next, scan backwards to trace out the upper boundary          
            for k=N:-1:1
                if validx(k)
                    % add new boundary point (avoid multiples)
                    d = Inf;
                    if numel(Bx)>0
                        d = sqrt((Xs(k)-Bx(end))^2+(Fs(k,iLow(k)+1)-By(end))^2);
                    end
                    if d>dtol
                        Bx = [Bx,Xs(k)];
                        By = [By,Fs(k,iLow(k)+1)];
                        Bsgn = [Bsgn,obj.ineq(Isrt(k,iLow(k)+1)).ineqSign.getDirection()];  
                        Bfcn = [Bfcn,Isrt(k,iLow(k)+1)];                        
                    end                    
                end
            end             
                        
            % finally, ensure that the boundary is closed
            Nb = numel(Bx);
            if Nb>2
                d = sqrt((Bx(1)-Bx(end))^2+(By(1)-By(end))^2);
                if d>dtol
                    Bx = [Bx,Bx(1)];
                    By = [By,By(1)];
                    Bsgn = [Bsgn,Bsgn(1)];
                    Bfcn = [Bfcn,Bfcn(1)];
                end
                obj.closed = true;
            end
            
            if feasrgn.feasrgn.DEBUG
                fprintf(1,'Number of boundary points: %d\n',Nb);
            end
            
            % removal of colinear segments
            P = [Bx; By];            
            Nb = numel(Bx);
            if Nb>2
                v = P(:,2:end)-P(:,1:end-1);
                v = v./(ones(2,1)*sqrt(sum(v.^2)));
                d = dot(v(:,1:end-1),v(:,2:end));
                idx = find(abs(abs(d)-1)<1e-11);
                if ~isempty(idx)
                    P(:,idx+1) = [];
                    Bfcn(idx+1) = [];
                end
            end
            obj.B = P;            
            Nb = size(P,2);         
            
            % generate inner points
%             generateInnerPts(obj, surfu, surfv)              

%{            
            fh = figure('Name','Feasible region demo');
            set(gca,'FontSize',10)
            % output constraints functions
            ph = plot(Xs,F);
            for i=1:nfcn
                set(ph(i),'LineStyle',obj.ineq(i).getLineStyle())
            end
            % output intersections S
            hold on
            plot(S(1,:),S(2,:),'o','MarkerSize',5,'Color',[0 0 0],'MarkerFaceColor',[1 1 1])
            hold off
            % output feasible region boundary
            hold on
            plot(obj.B(1,:),obj.B(2,:),'.-','LineWidth',2,'Color',[.75 0 0])
            hold off
            xmin = min(S(1,:)); xmax = max(S(1,:));
            ymin = min(S(2,:)); ymax = max(S(2,:));
            dx = (xmax-xmin)/10;
            dy = (ymax-ymin)/10;
            axis([xmin-dx/2 xmax+dx/2 ymin-dy/2 ymax+dy/2])
            grid on, box off
            xlabel(obj.ineq(1).iVarName)
            ylabel(obj.ineq(1).dVarName)
            title('Constraints functions & feasible region')
            legend(obj.legStr,'Location','eastoutside')
        
            % Output styled feasible region boundary
            %   the array Bfcn contains the indices of the inequalities
            
            segm = struct('hPlot',NaN,...
                'iFcn',NaN,...
                'iStart',NaN,...
                'iEnd',NaN,...
                'lineStyle','',...
                'lineColor',[]);            
            h = [find(diff(Bfcn)~=0),Nb];  
            segm(1).iFcn = Bfcn(1);
            segm(1).iStart = 1;
            segm(1).iEnd = h(1);
            segm(1).lineStyle = obj.ineq(Bfcn(1)).getLineStyle();
            segm(1).lineColor = get(ph(Bfcn(1)),'Color');            
            for k=1:length(h)-1
                segm(1+k).iFcn = Bfcn(h(k)+1);
                segm(1+k).iStart = h(k)+1;
                segm(1+k).iEnd = h(k+1);
                segm(1+k).lineStyle = obj.ineq(Bfcn(k)).getLineStyle();
                segm(1+k).lineColor = get(ph(Bfcn(k)),'Color');                
            end
            if feasrgn.feasrgn.DEBUG
                for i=1:length(segm)
                    segm(i)
                end
            end            
            figure, plot(obj.B(1,:),obj.B(2,:),'.-')
%}            
        end
            
        function set.queryPt(obj, q)      
            if isnumeric(q) && isreal(q) && numel(q)==2
                obj.queryPt = q(:);
            else
                error('Invalid query point')
            end
        end
        
        function res = isPtInside(obj, q)
            % returns true if query point q is inside the boundary B
            
%             figure
%             subplot(221),plot(obj.B(1,:),'.-')
%             subplot(223),plot(obj.B(2,:),'.-')
%             subplot(2,2,[2 4]),plot(obj.B(1,:),obj.B(2,:),'.-')     

            if ~obj.closed
                error('isPtInside requires a closed feasible region boundary')
            end
            if nargin>1
                if isnumeric(q) && isreal(q) && numel(q)==2
                    q = q(:);
                else
                    error('Input Q must be a 2-by-1 real array')
                end
                % set value for observer(s)
                obj.queryPt = q;                
            else
                if isempty(obj.queryPt)
                    error('Missing query point')
                else
                    q = obj.queryPt;
                end
            end              
                           
            Nb = size(obj.B,2);
            A = zeros(2*(Nb-1),2*(Nb-1));
            v = rand(2,1); % random ray
%             v = [1;1];
            for i=1:Nb-1
                A(1+2*(i-1) : 2*i, 1+2*(i-1) : 2*i) = [v,-(obj.B(:,i+1)-obj.B(:,i))];
            end            
            
%             det(A)
%             if abs(det(A))<1e-12
%                 error('Invalid boundary')
%             end                        

            % Tests if the query point q is inside of B
            % Create the coef. matrix b
            b = zeros(2*(Nb-1),1);
            for i=1:Nb-1
                b(1+2*(i-1) : 2*i) = obj.B(:,i)-q;
            end            
            % Solve for any intersection (all linear segments are tested
            % simultaneously)
            X = A\b;            
            % Count the number of intsects. betw. the ray v through q and B
            n = sum(X(1:2:end)>=0 & X(2:2:end)>=0 & X(2:2:end)<=1);
            % Odd number of intsects means that q is inside B
            res = (mod(n,2)~=0);            
        end
        
        function generateInnerPts(obj, surfu, surfv)           
            % surfu - number of segments in the u-direction
            % surfv - number of segments in the v-direction
            %
            % the u-direction is in the direction of the independent
            % variable
            
            u = linspace(0,1,surfu+1);
            v = linspace(0,1,surfv+1);  
            
            p = obj.p_;
            q = obj.q_;
                        
            t1 = feasrgn.feasrgn.parameterize(p);
            t2 = feasrgn.feasrgn.parameterize(q);
            
%             u = feasrgn.feasrgn.injectParams(u,t1);            
%             u = feasrgn.feasrgn.injectParams(u,t2);                        
            
            u = sort(u,'ascend');
            [idx1,h1] = feasrgn.feasrgn.mapParams(t1,u);
            [idx2,h2] = feasrgn.feasrgn.mapParams(t2,u);
            
            nu = length(u);
            nv = length(v);    
            % init outputs
            Gx_= zeros(nv,nu);
            Gy_ = zeros(size(Gx_));
            
            % generate inner points using Coon's bilinear interpolation
            % between the lower and upper boundary piecewise linear
            % functions
            P00 = p(:,1);
            P10 = p(:,end);
            P01 = q(:,1);
            P11 = q(:,end);
            
            for j=1:nv                
                % linear interp. in the v-direction
                P0v = p(:,1) + v(j)*(q(:,1)-p(:,1));
                P1v = p(:,end) + v(j)*(q(:,end)-p(:,end));                
                for i=1:nu           
                    % linear interp. in the u-direction
                    a = idx1(i);
                    b = idx2(i);
                    Pu0 = p(:,a)+(p(:,a+1)-p(:,a))*h1(i);
                    Pu1 = q(:,b)+(q(:,b+1)-q(:,b))*h2(i);  
                    % bilinear interpolation between points
                    Puv = (1-v(j))*Pu0+v(j)*Pu1+(1-u(i))*P0v+u(i)*P1v - ...
                        (1-u(i))*(1-v(j))*P00-v(j)*(1-u(i))*P01 - ...
                        u(i)*(1-v(j))*P10-u(i)*v(j)*P11;
                    Gx_(j,i) = Puv(1);
                    Gy_(j,i) = Puv(2);
                end
            end
            
            obj.Gx = Gx_;
            obj.Gy = Gy_;
%             figure
%             plot(p(1,:),p(2,:),'.-',q(1,:),q(2,:),'.-')                
%             hold on
%             plot(Gx,Gy,'.')
%             hold off
        end              
    end
    
    methods (Static)
        function [p] = removeColinear(p)
            % removes all colinear segments from the polyline P
            n = size(p,2);
            if n<3,return,end
            v = p(:,2:end)-p(:,1:end-1);
            v = v./(ones(2,1)*sqrt(sum(v.^2)));
            d = dot(v(:,1:end-1),v(:,2:end));
            idx = find(abs(abs(d)-1) < feasrgn.feasrgn.COLIN_TOL);
            if ~isempty(idx)
                p(:,idx+1) = [];
            end
        end
        function [t] = parameterize(p)
            % chordal parameterization of polyline P
            t = [0,cumsum(sqrt(sum((p(:,2:end)-p(:,1:end-1)).^2)))];
            t = t/t(end);
        end    
        function [u] = injectParams(u, t)
            % injects the polyline parameters 0<t<1 into the u-grid to
            % include the polyline definition points in the final surface.
            tol = 1e-11;
            nt = numel(t);
            for i=2:nt-1
                if ~any(abs(u-t(i)) < tol)
                    u = [u,t(i)];
                end
            end
        end
        function [idx,ulocal] = mapParams(t, u)
            % maps the parameter values u to boundary curve segment numbers 
            % to speed up later curve interpolation
            nt = length(t);
            nu = length(u);
            idx = zeros(1,nu);
            ulocal = zeros(1,nu);
            idx(1) = 1;
            ulocal(1) = 0;
            i = 2;
            for k=2:nu
                while u(k)>=t(i) && i<nt
                    i = i+1;
                end
                i = i-1;
                idx(k) = i;
                % compute local line parameter
                ulocal(k) = (u(k)-t(i))/(t(i+1)-t(i));
            end    
        end
    end   
end
