function [result,p_e,time,index_t]=ML_OA_MM(u,p,S,noise,SNR_Q,M,Q,d,p_t)
twpi = 2*pi;
epsilon=sqrt(2.2204e-10);%%smoothing factor
ux=u(1,:);
uy=u(2,:);
px=p(1);
py=p(2);
theta=atan2((ux-px),(uy-py));
J=size(S,2);
%% 
    sintheta1=sin(theta);
    costheta1=cos(theta);
    A = exp(-1i*twpi*d.'*[sintheta1;costheta1]);
    for i=1:Q
    S(i,:)=S(i,:)*sqrt(10^(SNR_Q(i)/10));
    end
    X=A*S+noise;
%%
    index_t=0; %%计数变量表示迭代次数
tic; %%开始计时
      Rxx=X*X'/J;%%求协方差矩阵
      p_tt=p_t+0.01;
      
while norm(p_t-p_tt)>=1e-4&&index_t<=5000%%相邻迭代点距离小于1e-8或者迭代次数大于500则停止迭代
          rt=zeros(1,Q); 
          for i=1:Q
                rt(i)=norm([u(:,i)-p_t;epsilon]);%%Smoothed distance
          end
          sintheta_t=(ux-p_t(1))./rt;
          costheta_t=(uy-p_t(2))./rt;
          At = exp(-1i*twpi*d.'*[sintheta_t;costheta_t]);
          
        Z2=(At'*At+epsilon*eye(Q))\(At'*Rxx*At)/(At'*At+epsilon*eye(Q)); 
        [~,eigenv]=eig(Z2);
        lambda_max_Z2=real(max(diag(eigenv)));
        Z2=(lambda_max_Z2*eye(Q)-Z2)*At';  

    
        Z1=(At'*At+epsilon*eye(Q))\At'*Rxx+Z2; 
        Amp_Z1=abs(Z1);
        real_Z1=real(Z1);
        imag_Z1=imag(Z1);
        Arg_Z1=atan2(imag_Z1,real_Z1);
        
        
     a_tt=0;
     b_tt=zeros(2,1);
    for q=1:Q
        A_11=zeros(3,3);
        b_11=zeros(3,1);
        for m=1:M
            [a,b]=cosine_surrogate(-2*Amp_Z1(q,m),Arg_Z1(q,m)-(twpi*d(:,m).'*[sintheta_t(q);costheta_t(q)]));
            A_11=blkdiag(a*(twpi)^2*d(:,m)*d(:,m).',0)+A_11;
            b_11=-2*a*Arg_Z1(q,m)*twpi*[d(:,m) ;0]-b*twpi*[d(:,m) ;0]+b_11;
        end
             [a_t,b_t,~]=Rayleigh_Quotient_Surrogate(A_11,b_11,[u(:,q)-p_t;epsilon],0);
             a_tt=a_t+a_tt;
             b_tt=-b_t(1:2)-2*a_t*u(:,q)+b_tt;
             
           
    end
%%    
p_tt=p_t;    
% p_t=-b_tt/(2*a_tt);
p_t = nearestPointOnSquare(-b_tt/(2*a_tt));


%%
      index_t=index_t+1;%%迭代加一
        cc(:,index_t)=p_tt;
     ob_value(index_t)=objective_ML_OA(Rxx,u,ux,uy,p_t,d,M,twpi,Q,epsilon);
end
  norm(find(ob_value(1:end-1)-ob_value(2:end)<=0));  %%判断目标函数是否下降，用于验证代码正确性,在迭代次数较大时，由于精度会导致不成立  
%     semilogy(ob_value,'-k')  
   p_e=p_t;   
time=toc;
result=sum((p_e-p).^2);