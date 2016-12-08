%{
Finds the piecewise linear feasible region from a set of constraints 
functions (inequalities) derived within the actual problem domain.

    y 'sign' f1(x)
    y 'sign' f2(x)
    :
    y 'sign' fn(x)

where 'sign' can be any inequality sign '<'|'<='|'>'|'>=' and f(x) is the 
constraints function which is either constant, linear or possibly nonlinear

Description
The algorithm makes no assumptions on the specified constraints functions 
which can be either constant, linear or possibly nonlinear. The constraints 
functions are initially represented by anonymous functions which will be 
converted into piecewise linear functions within the specified interval 
xmin<=x<=xmax for the independent variable.
The resolution of the discrete grid is controlled by (dx).

It works out the feasible region by first sorting the sampled functions 
followed by detection of their intersections which will be injected into 
the discrete x-grid so that they can be included in the final feasible 
region boundary. The final feasible region is then found by scanning along 
the discrete x-grid for valid (non-conflicting) constraints.

The output is ultimately a closed piecewise linear boundary B representing 
the feasible region to be used in solving 2-dimensional optimization 
problems.

Possibly multiple feasible regions is not handled which could yield
unexpected results.

Dependencies:
Package: feasrgn

Version : 1.2
Date    : 12/2016
Author  : Jari Repo, University West, jari.repo@hv.se
--------------------------------------------------------------------------
%}
clear all, close all, clc
format short

% ------------------------------------------------------------------------
% Include package feasrgn
% ------------------------------------------------------------------------
import feasrgn.*

% ------------------------------------------------------------------------
% Define the constraints functions (stored in cell-array cfcn)
% ------------------------------------------------------------------------
Ex = 2;     % select example 1-4 
[cfcn,opts] = LoadFeasRgnExample( Ex );

% ------------------------------------------------------------------------
% Generate the piecewise linear feasible region
% ------------------------------------------------------------------------
% To inspect individual inequality objects: frgn.ineq( index )
frgn = feasrgn(cfcn,opts);

if ~frgn.closed
    error('No closed feasible region')
end

% ------------------------------------------------------------------------
% Init figure and plot the feasible region
% ------------------------------------------------------------------------
fontSize = 9;

fh = figure('Name','Piecewise linear feasible region (v1.2)');
set(gca,'FontSize',fontSize)
title('Constraints functions & feasible region')
box off,grid off
xlabel(frgn.ineq(1).iVarName)
ylabel(frgn.ineq(1).dVarName)

obj = feasrgnplot(frgn,...
    'LineStyle','none',...
    'LineWidth',2,...
    'FillColor',[0 .45 .85],...  
    'FaceAlpha',.9,...
    'NodeLabel','S',...
    'NodeShape','o',...
    'NodeSize',5,...
    'NodeColor',[0 0 .75],...
    'NodeFontSize',8,...
    'DisplayNodes',true,...  
    'DisplayNodeValues',true,...
    'DisplayQueryPoint',false,...
    'NodeValueFormat','(%.2f; %.2f)',...
    'DisplayConstraints',true,...
    'DisplayLegend',true,...
    'DisplayInnerPoints',true,...
    'UPointCount',30,...
    'VPointCount',20);

disp('Any key...')
pause

% ------------------------------------------------------------------------
% Modify feasible region plot properties
% ------------------------------------------------------------------------
% also possible to dynamically change any of the above properties

obj.nodelabel = 'P';
obj.nodeshape = 's';    % square
obj.nodecolor = 'green';
obj.nodevalueformat = '[%.2f, %.2f]';

disp('Any key...')
pause

% ------------------------------------------------------------------------
% Point inside boundary test
% ------------------------------------------------------------------------
% Test if query point q=[x;y] is inside the feasible region boundary B
obj.displayquerypoint = true;
q = [-.5;2];
res = frgn.isPtInside(q);
if res==true
    fprintf(1,'Point (%f; %f) is inside the boundary\n',q');
else
    fprintf(1,'Point (%f; %f) is outside the boundary\n',q');
end            

pause(1)

frgn.queryPt = [0;2.5];   % this will change the query point marker on-screen
res = frgn.isPtInside();

q = frgn.queryPt;
if res==true
    fprintf(1,'Point (%f; %f) is inside the boundary\n',q');
else
    fprintf(1,'Point (%f; %f) is outside the boundary\n',q');
end            

% ------------------------------------------------------------------------
% Generate inner boundary points
% ------------------------------------------------------------------------
% The generated point set can be used to evaluate the objective function
% for the current optimization problem at hand

%   frgn.generateInnerPts(30,40)
%   frgn.Gx
%   frgn.Gy


% % Output the text '+feasrgn' for the File Exchange cover image
% hold on
% text(mean(frgn.Gx(:)),mean(frgn.Gy(:)),1,'+feasrgn',...
%     'Clipping','on','FontSize',14,'Color',[1 1 1],'FontName','Cambria',...
%     'FontWeight','normal','HorizontalAlignment','center',...
%     'VerticalAlignment','bottom','Visible','on')
% hold off
% PNG_Export(gcf,20,15,'feasible_region_1_2');
