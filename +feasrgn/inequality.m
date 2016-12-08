classdef inequality < handle
% inequality Represents an inequality of the form y 'sign' f(x)
% 
% where 'sign' is one of the following
%   '<'  : less than
%   '<=' : less than or equal to
%   '>'  : greater than
%   '>=' : greater than or equal to
%
% The function f(x) is given as an anonymous function handle which can
% define either a constant, linear or nonlinear function.
% 
% USAGE: ineq = inequality(number, dVarName, ineqSignStr, hFcn) 
%
%   number - index used in the right-hand-side name
%   dVarName - name of dependent variable, for example 'y'
%   ineqSignStr - inequality sign symbol, for example '<','<=','>' or '>='
%   hFcn - anonymous single-variable function handle, for example:
%               @(x) x.^2
%               @(x) 2*ones(size(x)
%               @(t) sin(pi*t/2).^2
% 
% inequality Properties:
%   iVarName (read-only)
%   dVarName (read-only)
%   number (read-only)
%   ineqSign (read-only)
%   hFcn (read-only)
%   lhsName
%   rhsName
% 
% inequality Methods:
%   setDefaultNames
%   getSignature
%   getLineStyle
%   isNeg
%   isPos
% 
% Dependencies:
%   feasrgn.ineqsign
% 
% Features:
%   -the name of the independent variable is obtained directly from the 
%    specified anonymous function
%   -possible to override the predefined left-hand-side and right-hand-side
%   signatures by using the corresponding property-set methods.

%{
    --- Change log ---

    Version 1.2, JRE
    -renamed class to 'inequality'
    -removed the auto-indexing
    -replaced lt() and gt() with isNeg() and isPos() respectively
    -the signature of the inequality is obtained by getSignature()
    -default lhsName and rhsName are set by setDefaultNames()
    -the logics related to the inequality sign is now handled by the 
    enumeration class ineqsign

    Version 1.1, JRE
    -New methods:
    lt() : return true if the inequality sign is '<' or '<='
    gt() : return true if the inequality sign is '>' or '>='
    
    Version : 1.2
    Date    : 2016-12-06
    Author  : Jari Repo, University West, jari.repo@hv.se
%}

    properties
        lhsName = '';       % left-hand-side signature
        rhsName = '';       % right-hand-side signature       
    end
    properties (SetAccess = private)
        number = NaN;                       % inequality index 1,2,...
        ineqSign = feasrgn.ineqsign.empty;  % inequality sign
        hFcn = function_handle.empty;       % function handle
        iVarName = '';      % name of independent variable, 'x'
        dVarName = '';      % name of dependent variable, 'y'  
    end
    properties (Constant, Hidden = true)
        solidLine = '-';            % line style for '<=' and '>='
        dashedLine = '--';          % line style for '<' and '>'
        fsymbol = 'f';              % right-hand-side function label
    end
    methods
        function obj = inequality(number, dVarName, ineqSignStr, hFcn)
            if nargin>3
                if isnumeric(number) && isreal(number) && ~isempty(number) && ...
                        numel(number)==1 && number>=0 && number==fix(number)
                    obj.number = number;
                else
                    error('Input NUMBER must be a real number')
                end
                if isa(dVarName,'char') && ~isempty(strtrim(dVarName))
                    obj.dVarName = strtrim(dVarName);
                else
                    error('Input DVARNAME must be of type ''char''')
                end
                if isa(ineqSignStr,'char')
                    obj.ineqSign = feasrgn.ineqsign.enum(ineqSignStr);
                else
                    error('Input INEQSIGNSTR must be of type ''char''')
                end
                if isa(hFcn,'function_handle')
                    % get the name of the independent variable
                    % https://se.mathworks.com/matlabcentral/answers/121605-regexp-parentheses-string-parsing-issue
                    % obj.iVarName = regexp(func2str(hFcn),'(?<=\()\S','match','once');
                    obj.hFcn = hFcn;
                    obj.iVarName = strtrim(regexp(func2str(obj.hFcn),...
                        '(?<=\().+?(?=\))','match','once'));
                    obj.setDefaultNames();
                else
                    error('Input HFCN must be a function handle')
                end
            else
                error('Missing inputs')
            end
        end
        function setDefaultNames(obj)
            % sets the default lhsName and rhsName
            obj.lhsName = obj.dVarName;
            obj.rhsName = sprintf('%s_{%d}(%s)',...
                feasrgn.inequality.fsymbol,...
                obj.number,...
                obj.iVarName);
        end
        function value = getSignature(obj)
            % returns the inequality signature, for example 'y < f_{1}(x)' 
            % to be shown in figure legends
            value = sprintf('%s %s %s',...
                obj.lhsName,obj.ineqSign.getSymbol(),obj.rhsName);
        end
        function value = getLineStyle(obj)
            % returns the appropriate line style when plotting the 
            % constraints function
            if obj.ineqSign==feasrgn.ineqsign.LT || ...
                    obj.ineqSign==feasrgn.ineqsign.GT
                value = feasrgn.inequality.dashedLine;
            else
                value = feasrgn.inequality.solidLine;
            end
        end
        function value = isNeg(obj)
            % returns true if inequality sign is '<' or '<='
            value = ( obj.ineqSign.getDirection()==-1 );
        end
        function value = isPos(obj)
            % returns true if inequality sign is '>' or '>='
            value = ( obj.ineqSign.getDirection()==1 );
        end
        function set.lhsName(obj, value)
            % sets the left-hand-side signature
            if isa(value,'char') && ~isempty(value)
                obj.lhsName = value;
            else
                error('Length of LHSNAME must be > 0')
            end
        end        
        function set.rhsName(obj, value)
            % sets the right-hand-side signature
            if isa(value,'char') && ~isempty(value)
                obj.rhsName = value;
            else
                error('Length of RHSNAME must be > 0')
            end
        end      
    end
end
