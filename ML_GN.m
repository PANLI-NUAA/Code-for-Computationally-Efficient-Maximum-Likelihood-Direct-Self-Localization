function [result,p_e,time,numIter]=ML_GN(u,p,S,noise,SNR_Q,M,Q,d,p_t)
twpi = 2*pi;
ux=u(1,:);
uy=u(2,:);
px=p(1);
py=p(2);
theta=atan2((ux-px),(uy-py));
J=size(S,2);
tic;
%% 
    sintheta1=sin(theta);
    costheta1=cos(theta);
    A = exp(-1i*twpi*d.'*[sintheta1;costheta1]);
%%    
    
    for i=1:Q
    S(i,:)=S(i,:)*sqrt(10^(SNR_Q(i)/10));
    end
    X=A*S+noise;
%%
        Rxx=X*X'/J;%%求协方差矩阵  


 fun = @(x) myFunction(x,Rxx,u,M,d);

% 初始点
x0 = p_t;

% 配置选项，选择 BFGS 算法
options = optimoptions('fminunc', ...
    'FunctionTolerance', 1e-12,  ...
    'StepTolerance', 1e-4,     ...  
    'MaxIterations',5000,    ...   
    'Display', 'off');      
   


% 调用 fminunc
[p_e, fval,exitflag,output] = fminunc(fun, x0, options);
numIter = output.iterations;
% disp(['最优解: ', mat2str(x)]);
% disp(['目标函数值: ', num2str(fval)]);

      
  time=toc;    
% p_e = nearestPointOnSquare(p_e);
result=sum((p_e-p).^2);