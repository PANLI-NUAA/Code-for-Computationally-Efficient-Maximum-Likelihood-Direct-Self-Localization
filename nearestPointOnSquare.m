function q = nearestPointOnSquare(p)
    % p: 1x2 向量 [x y]
    % q: 最近的点在正方形边界上的坐标
    
    xmin = -2000; xmax = 2000;
    ymin = -2000; ymax = 2000;
    x = p(1); y = p(2);

    % 限制点到正方形范围
    x_clamped = min(max(x, xmin), xmax);
    y_clamped = min(max(y, ymin), ymax);

    % 如果点在外部，最近点就是 clamp 后的点
    if x < xmin || x > xmax || y < ymin || y > ymax
        q = [x_clamped, y_clamped].';
        return;
    end
    q=[x,y].';
    % 如果点在正方形内部，投影到最近的边
%     dx_left   = abs(x - xmin);
%     dx_right  = abs(x - xmax);
%     dy_bottom = abs(y - ymin);
%     dy_top    = abs(y - ymax);
% 
%     [~, idx] = min([dx_left, dx_right, dy_bottom, dy_top]);
% 
%     switch idx
%         case 1, q = [xmin, y].';
%         case 2, q = [xmax, y].';
%         case 3, q = [x, ymin].';
%         case 4, q = [x, ymax].';
%     end
end
