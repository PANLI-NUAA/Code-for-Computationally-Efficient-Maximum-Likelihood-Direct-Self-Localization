function fval = myFunction(x,R,u,M,d)
    % 输入：x 是变量向量
    % 输出：fval 是目标函数值

    % 复杂目标函数定义
    theta=atan2((u(1,:)-x(1)),(u(2,:)-x(2)));

   sintheta1=sin(theta);
    costheta1=cos(theta);
    A = exp(-1i*2*pi*d.'*[sintheta1;costheta1]);
    
    fval =real(trace((eye(M)-A/(A'*A)*A')*R));
end
