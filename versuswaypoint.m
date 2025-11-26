clear all
close all
maxNumCompThreads('automatic');
NsePwrVecdB =3:8; 
J=126;
Mc=1000;
u=[2000,4000; 6000 2000; -3000 4150; 2000 -3500; -1725 -5600; -3500 2000; -1425 3000 ; -3000 1775].';
Q=length(u(1,:));

M2=10; %%阵列的阵元数目

randn('state',1);
S_r=randn(Q,J);
S_i=randn(Q,J);
S=(S_r+1i*S_i)/sqrt(2);
% mean(mean(abs(S).^2))%%该代码用于验证平均功率是否为1
for n=1:Mc
N_r{n}=randn(M2,J);
N_i{n}=randn(M2,J);
noise{n}=(N_r{n}+1i*N_i{n})/sqrt(2);
%  mean(mean(abs(noise{n}).^2)); %%该代码用于验证平均功率是否为1
end
p=[34.43,53.875].';
dd=0.5;
d=0:dd:(M2-1)*dd;

d(2,:)=zeros(1,M2);    
% for i=1:M2
%     d(1,i)=cos(2*pi/M2*i);
%     d(2,i)=sin(2*pi/M2*i);
% end
 SNR=-20;
Number_of_anchors=NsePwrVecdB(end);
        r0=zeros(Number_of_anchors,1); 
        for i=1:Number_of_anchors
                r0(i)=norm(u(:,i)-p);%%求真实距离
        end       
        kappa = 10^(SNR./10)*Number_of_anchors/sum(r0.^(-2));
        SNR_Q1 = kappa*r0.^(-2);
        SNR_Q1=10*log10(SNR_Q1);%%各个辐射源的SNR
rand('state',4);
for j=1:length(NsePwrVecdB)
    SNR=-20;%%平均SNR，接下来求各个anchor辐射源的SNR。假设SNR随着1/距离的平方衰减
%     SNR=10;

Number_of_anchors=NsePwrVecdB(j);
SNR_Q=SNR_Q1(1:Number_of_anchors);

        %%10*log10(mean(10.^(SNR_Q/10))) %%验证平均SNR是否正确
parfor i=1:Mc 
range=2000;
S2=S(1:Number_of_anchors,1:J);
noise2=noise{i};
    [result_WSSF(i),p_e_WSSF(:,i),time_WSSF(i)]=WSSF(u(:,1:Number_of_anchors),p,-range,range,-range,range,10,S2,noise2,SNR_Q,M2,d,Number_of_anchors);
    [result_ML_OA(i),p_e_ML_OA(:,i),time_ML_OA(i)]=ML_OA(u(:,1:Number_of_anchors),p,-range,range,-range,range,10,S2,noise2,SNR_Q,M2,d,Number_of_anchors);   
     [~,p_e2,time_zz]=ML_OA(u(:,1:Number_of_anchors),p,-range,range,-range,range,100,S2,noise2,SNR_Q,M2,d,Number_of_anchors);    
       [result_ML_OA_MM2(i),p_PLE_ML_OA_MM2(:,i),time_ML_OA_MM12,ttt]=ML_OA_MM(u(:,1:Number_of_anchors),p,S2,noise2,SNR_Q,M2,Number_of_anchors,d,p_e2);
        [result_ML_GN2(i),p_PLE_ML_GN2(:,i),time_ML_GN2,~]=ML_GN(u(:,1:Number_of_anchors),p,S2,noise2,SNR_Q,M2,Number_of_anchors,d,p_e2);
    time_ML_OA_MM2(i)=time_ML_OA_MM12+time_zz;
     time_GN_OA2(i)=time_ML_GN2+time_zz;  
end
f_WSSF(j)=sqrt(mean(result_WSSF));
f_ML_OA(j)=sqrt(mean(result_ML_OA));
f_ML_OA_MM2(j)=sqrt(mean(result_ML_OA_MM2));
f_ML_GN2(j)=sqrt(mean(result_ML_GN2));

time_WSSF1(j)=mean(time_WSSF);
time_ML_OA1(j)=mean(time_ML_OA);
time_ML_OA_MM122(j)=mean(time_ML_OA_MM2);
time_GN_OA122(j)=mean(time_GN_OA2);
end

for j=1:length(NsePwrVecdB)  
    %%平均SNR，接下来求各个anchor辐射源的SNR。假设SNR随着1/距离的平方衰减
%     SNR=10;
Number_of_anchors=NsePwrVecdB(j);
SNR_Q=SNR_Q1(1:Number_of_anchors);

     CRLB(j)=(SPEB(M2,SNR_Q(1:Number_of_anchors),J,p,u(:,1:Number_of_anchors),d,Number_of_anchors));

end



figure(1)
plot(NsePwrVecdB,f_WSSF, 's', 'MarkerEdgeColor', '#77AC30', 'MarkerFaceColor', 'none', 'MarkerSize', 8, 'LineWidth', 1.2)
hold on
plot(NsePwrVecdB,f_ML_OA, 'o', 'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor', 'none', 'MarkerSize', 8, 'LineWidth', 1.2)
hold on
plot(NsePwrVecdB,f_ML_GN2,'*','Color', '#0000FF', 'MarkerEdgeColor', '#0000FF', 'MarkerFaceColor', '#0000FF', 'MarkerSize', 8, 'LineWidth', 1.2)
hold on
plot(NsePwrVecdB,f_ML_OA_MM2,'x', 'MarkerEdgeColor', '#A2142F', 'MarkerFaceColor', 'none', 'MarkerSize', 8,'LineWidth', 1.2)
hold on

plot(NsePwrVecdB,CRLB,'-k','linewidth',0.5)
legend('WSSF','ML (exhaustive search)','ML-QN','ML-MM (proposed)','CRLB','Fontsize',12)
xlabel('SNR [dB]','Fontsize',12)
ylabel('RMSE [m]','Fontsize',12)

