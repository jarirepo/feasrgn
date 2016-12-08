function [cfcn,opts] = LoadFeasRgnExample(n)
% LoadFeasRgnExample Predefined feasible region examples
% 
% Inputs
%   n - example number (1-4)
% 
% Outputs
%   cfcn - a cell array holding all constraints functions (inequalities)
%   opts - options structure with fields 'xmin','xmax','dx'

%{
    Version : 1.0
    Date    : 2016-12-05
    Author  : Jari Repo, University West, jari.repo@hv.se
%}
    switch n
        case 1 % Example #1
            cfcn = { 'y','<',@(x).675-x/2; ...
                     'y','>',@(x) 2*x.*(x-1).^2; ...
                     'y','<',@(x) sin(pi/2*x).^2; ...
                     'y','<=',@(x) cos(pi/2*x).^2 };
                 
            % Specify the x-limits [xmin,xmax] and the resulution (dx)
            opts = struct('xmin',.25,'xmax',.75,'dx',0.01);
            
            % It is possible to specify an empty options structure.
            % In that case, the feasrgn class will then find suitable 
            % values for xmin,xmax,dx            
            % opts = {};
            
        case 2 % Example #2
            cfcn = { 'y','>', @(x) .65*(x-1).^2-3; ...
                     'y','<=',@(x) -1.25*x.^2+6; ...
                     'y','>=', @(x) .8*(x+1.25).^3-1.025*(x+.5); ...
                     'y','>', @(x) .5*ones(size(x)); ...
                     'y','<=',@(x) 3.5*ones(size(x)); ...
                     'y','<=',@(x) -3*x+2; ...
                     'y','<=',@(x) 3-.5*x };
                
            opts = struct('xmin',-3,'xmax',3,'dx',0.05);
%             opts = {};
            
        case 3 % Example #3
            % Define 3 constraints functions:
            %
            %   r > 0
            %   r < -4t(t-1/2)
            %   r > t
            %
            % Here we use different names for the independent and dependent
            % variables, 't' and 'r' respectively
            
            cfcn = { 'r','>', @(t) zeros(size(t)); ...
                     'r','<', @(t) -4*t.*(t-1/2); ...
                     'r','>', @(t) t };
                        
            opts = struct('xmin',-1,'xmax',1,'dx',.01);
            
        case 4 % Example #4
            % Linear inequalities of type f(x)=a0+a1*x       
            cfcn = { 'y','>', @(x) .1-x/10; ...
                     'y','<', @(x) 2*x; ...
                     'y','<', @(x) 1.5*x; ...
                     'y','<=',@(x) (1/2)*(1-x); ...
                     'y','<', @(x) x;
                     'y','<', @(x) 1-x; ...
                     'y','>=',@(x) .1+x/10 };
                
            opts = struct('xmin',-1,'xmax',2,'dx',.01);
            
        otherwise
            error('No inequality set defined')
    end
end
