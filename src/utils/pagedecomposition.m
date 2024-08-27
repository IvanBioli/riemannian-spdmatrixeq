classdef pagedecomposition
    % pagedecomposition - Decomposition of a three-dimensional array
    properties
        r
        decomposition
    end
    methods
        function obj = pagedecomposition(G, varargin)
            % pagedecomposition - Class object definition
            
            r = size(G, 3); obj.r = r;
            obj.decomposition = cell(r, 1);
            for i = 1:r
                obj.decomposition{i} = decomposition(G(:, :, i), varargin{:});
            end
        end
        function x = pagemldivide(obj, h)
            % pagemldivide - mldivide
            
            x = zeros(size(h));
            for i = 1:obj.r
                x(:, :, i) = obj.decomposition{i} \ h(:, :, i);
            end
        end
    end
end