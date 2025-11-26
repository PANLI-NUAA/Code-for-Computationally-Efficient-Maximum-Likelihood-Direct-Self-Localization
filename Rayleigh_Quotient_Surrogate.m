function [a_t,b_t,const]=Rayleigh_Quotient_Surrogate(A,b,x_t,mode)
%%The job of this function is to calculate the coefficients in Theorem 2.
%%Input A, b, and x_t correspond to A, b, and x_t in Theorem 2.
%%Output a_t and b_t correspond to a_t and b_t in Theorem 2.
%%If mode=1, we calculate the const. If mode=0, we did not compute the const.

B=[A,0.5*b;0.5*b.',0];

D=length(x_t); %% Calculate the dimension. For the three-dimensional case, D=3. For the two-dimensional case, D=2.  

lambda_max=norm(B);%% Calculate the spectrum norm (maximum eigenvalue).

lambda=2*lambda_max; %% In our simulations, $\epsilon$ is set to $\lambda_{max}(\mathbf B)$.

C=lambda*eye(D+1)-B; %% Compute matrix C;
invC=inv(C);

P=blkdiag([norm(x_t)*eye(D)],1); 
invP=blkdiag([1/norm(x_t)*eye(D)],1);

y_t=invP*C*invP*[x_t;1];%% Compute y_t;

widetildey_t=y_t(1:D);
consy_t=y_t(end);

widetildeD=invC(1:D,1:D); %%Compute $\widetilde{\mathbf D}$.
d=invC(1:D,end);
d_const=invC(end,end);

Z=consy_t*d.'*widetildey_t;

if Z>=0
a_t=Z/norm(x_t)+widetildey_t.'*widetildeD*widetildey_t;
b_t=-2*widetildey_t;
else
    a_t=widetildey_t.'*widetildeD*widetildey_t;
    b_t=-2*widetildey_t+2*Z/norm(x_t)*x_t;
end

if mode==1%%Calculate the nonsense constant.
    if Z>=0
const=d_const*consy_t^2+2*lambda-2*consy_t+norm(x_t)*Z;
    else 
        const=d_const*consy_t^2+2*lambda-2*consy_t;
    end
else
    const=0;
end


