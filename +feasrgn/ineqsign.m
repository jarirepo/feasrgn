classdef ineqsign
% ineqsign  Inequality sign enumeration
% 
% ineqsign Methods:
%   getDirection
%   getSymbol
%   ineqsign.enum

%{
    Version : 1.0
    Date    : 2016-12-05
    Author  : Jari Repo, University West, jari.repo@hv.se
%}
    enumeration
        LT, LE, GT, GE
    end      
    methods
        function value = getDirection(obj)            
            if obj.lt() || obj.le()
                value = -1;
            else
                value = 1;
            end            
        end
        function value = getSymbol(obj)
            switch obj
                case feasrgn.ineqsign.LT, value = '<';
                case feasrgn.ineqsign.LE, value = '\leq';
                case feasrgn.ineqsign.GT, value = '>';
                case feasrgn.ineqsign.GE, value = '\geq';
            end
        end
        function value = lt(obj)
            value = ( obj == feasrgn.ineqsign.LT );
        end
        function value = le(obj)
            value = ( obj == feasrgn.ineqsign.LE );
        end
        function value = gt(obj)
            value = ( obj == feasrgn.ineqsign.GT );
        end
        function value = ge(obj)
            value = ( obj == feasrgn.ineqsign.GE );
        end        
    end
    methods (Static)
        function obj = enum(str)
            switch lower(str)
                case {'<','lt'}
                    obj = feasrgn.ineqsign.LT;
                case {'<=','=<','le','leq','\leq'}
                    obj = feasrgn.ineqsign.LE;
                case {'>','gt'}
                    obj = feasrgn.ineqsign.GT;
                case {'>=','=>','ge','geq','\geq'}
                    obj = feasrgn.ineqsign.GE;
                otherwise
                    error('Inequality sign ''%s'' is not recognized')
            end
        end
    end
end
