function [result,p_e,time]=ML_OA(u,p,xmin,xmax,ymin,ymax,acc,S,noise,SNR_Q,M,d,Q)
radeg = 180/pi;
twpi = 2*pi;
xx=xmin:acc:xmax;
yy=ymin:acc:ymax;
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
    tic;
      Rxx1=X*X'/J;
      lx=length(xx);
ly=length(yy);
temp1=zeros(lx,ly);
    for i=1:lx
        for j=1:ly  
          for q=1:Q
                rt(q)=norm([u(:,q)-[xx(i);yy(j)]]);%%Smoothed distance
          end
          sintheta=(ux-xx(i))./rt;
          costheta=(uy-yy(j))./rt;
            A=exp(-1i*2*pi*d.'*[sintheta ;costheta]);            
temp1(i,j)=-trace((A*pinv(A))*Rxx1);
        end
    end
SP1=real(temp1);

[x,y]=find(SP1==min(min(SP1)));
x=(x-1)*acc+xmin;
y=(y-1)*acc+ymin;
p_e=[x y].';
% h=mesh(xx,yy,1./SP1');
% set(h,'Linewidth',2);
% xlabel('x/m');
% ylabel('y/m');
% zlabel('幅度/dB');
% view(2)
% axis([xmin xmax ymin ymax 0 1]);
time=toc;
result=sum((p_e-p).^2);