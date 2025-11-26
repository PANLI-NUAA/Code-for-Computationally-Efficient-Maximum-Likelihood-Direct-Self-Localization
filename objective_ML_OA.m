function ob_value=objective_ML_OA(Rxx,u,ux,uy,p_t,d,M,twpi,Q,epsilon)
          rt=zeros(1,Q); 
          for i=1:Q
                rt(i)=norm([u(:,i)-p_t;epsilon]);%%Smoothed distance
          end
          sintheta_t=(ux-p_t(1))./rt;
          costheta_t=(uy-p_t(2))./rt;
          A = exp(-1i*twpi*d.'*[sintheta_t;costheta_t]);
      ob_value=real(trace((eye(M)-A/(A'*A+epsilon*eye(Q))*A')*Rxx));