%$%
% This method should be inside GraphSLAM class.
% NOTE: should return a (x, 1) matrix

function idx = find_iMat_idx(obj, num, str)
    % get indices of the info matrix
    % num: can be a single num or a list
    % str: 'x': pose; 'm': maps
    
    % GLOBAL VARIABLES
    numX = size(obj.x, 2);
    
    if str == 'x'
        len = size(num(:), 1);
        if len == 1
            idx = [3 * num - 2, 3 * num - 1, 3 * num];
            return
        end
        
        list = reshape(num, [len, 1]);
        c1 = 3 * list - 2;
        c2 = c1 + 1;
        c3 = c2 + 1;
        c = [c1, c2, c3];
        c1 = c';
        idx = c1(:);
        
    end
    
    if str == 'm'
        bias = 3 * numX;
        len = size(num(:), 1);
        if len == 1
            idx = [3 * num - 2, 3 * num - 1, 3 * num];
            idx = idx + bias;
            return
        end
        
        list = reshape(num, [len, 1]);
        c1 = 3 * list - 2;
        c2 = c1 + 1;
        c3 = c2 + 1;
        c = [c1, c2, c3];
        c1 = c';
        idx = c1(:);
        idx = idx + bias;
    end
    
end