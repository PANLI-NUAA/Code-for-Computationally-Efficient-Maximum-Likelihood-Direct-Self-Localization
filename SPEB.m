function f=SPEB(M,SNR,J,p,u,d,Q)
twpi = 2*pi;
ux=u(1,:);
uy=u(2,:);
px=p(1);
py=p(2);
theta=atan2((ux-px),(uy-py));
    sintheta1=sin(theta);
    costheta1=cos(theta);
    A = exp(-1i*twpi*d.'*[sintheta1;costheta1]);    
    P_perp_A=eye(M)-A*((A'*A)\A');

    for i=1:Q
        for j=1:M
            D(j,i)=(-1i*twpi*d(:,j).'*[costheta1(i);-sintheta1(i)])*exp(-1i*twpi*d(:,j).'*[sintheta1(i);costheta1(i)]);
        end
    end
    
   CRB1= real(2*J*(D'*P_perp_A*D).*diag(10.^(SNR./10)));
  
   for i=1:Q
       r=norm(p-u(:,i));
   T(:,i)=[costheta1(i) -sintheta1(i)].'./r;
   end
   
CRB2=inv(T*CRB1*T.');


f=sqrt(trace(CRB2));

