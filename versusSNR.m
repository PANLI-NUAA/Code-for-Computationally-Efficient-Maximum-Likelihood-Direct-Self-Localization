clear all
close all
maxNumCompThreads('automatic');
NsePwrVecdB =-20:2.5:20; 
J=128;
Mc=1000;
u=[2000,4000; 6000 2000; -3000 4150].';
Q=length(u(1,:));

M2=5; %%阵列的阵元数目

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

rand('state',4);
p_e2r=400*rand(2,Mc)-200;
for j=1:length(NsePwrVecdB)
    SNR=NsePwrVecdB(j);%%平均SNR，接下来求各个anchor辐射源的SNR。假设SNR随着1/距离的平方衰减
%     SNR=10;
        r0=zeros(Q,1); 
        for i=1:Q
                r0(i)=norm(u(:,i)-p);%%求真实距离
        end       
        kappa = 10^(SNR./10)*Q/sum(r0.^(-2));
        SNR_Q = kappa*r0.^(-2);
        SNR_Q=10*log10(SNR_Q);%%各个辐射源的SNR
        %%10*log10(mean(10.^(SNR_Q/10))) %%验证平均SNR是否正确
parfor i=1:Mc 
range=2000;
    [result_ML_OA(i),p_e_ML_OA(:,i),time_ML_OA(i)]=ML_OA(u,p,-range,range,-range,range,10,S,noise{i},SNR_Q,M2,d,Q);   
    [result_WSSF(i),p_e_WSSF(:,i),time_WSSF(i)]=WSSF(u,p,-range,range,-range,range,10,S,noise{i},SNR_Q,M2,d,Q);
     [~,p_e2,time_zz]=ML_OA(u,p,-range,range,-range,range,100,S,noise{i},SNR_Q,M2,d,Q);    
       [result_ML_OA_MM2(i),p_PLE_ML_OA_MM2(:,i),time_ML_OA_MM12,index_OA_MM(i)]=ML_OA_MM(u,p,S,noise{i},SNR_Q,M2,Q,d,p_e2);
        [result_ML_GN2(i),p_PLE_ML_GN2(:,i),time_ML_GN2,index_ML_GN(i)]=ML_GN(u,p,S,noise{i},SNR_Q,M2,Q,d,p_e2);
    time_ML_OA_MM2(i)=time_ML_OA_MM12+time_zz;
     time_GN_OA2(i)=time_ML_GN2+time_zz;  
end
iteration_OA_MM(j)=mean(index_OA_MM);
iteration_ML_GN(j)=mean(index_ML_GN);
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
     SNR=NsePwrVecdB(j);
        r0=zeros(Q,1); 
        for i=1:Q
                r0(i)=norm(u(:,i)-p);%%求真实距离
        end       
        kappa = 10^(SNR./10)*Q/sum(r0.^(-2));
        SNR_Q = kappa*r0.^(-2);
        SNR_Q=10*log10(SNR_Q);%%各个辐射源的SNR

     CRLB(j)=(SPEB(M2,SNR_Q,J,p,u,d,Q));

end
figure(1)
semilogy(NsePwrVecdB,f_WSSF, 's', 'MarkerEdgeColor', '#77AC30', 'MarkerFaceColor', 'none', 'MarkerSize', 8, 'LineWidth', 1.2)
hold on
semilogy(NsePwrVecdB,f_ML_OA, 'o', 'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor', 'none', 'MarkerSize', 8, 'LineWidth', 1.2)
hold on
semilogy(NsePwrVecdB,f_ML_GN2,'*','Color', '#0000FF', 'MarkerEdgeColor', '#0000FF', 'MarkerFaceColor', '#0000FF', 'MarkerSize', 8, 'LineWidth', 1.2)
hold on
semilogy(NsePwrVecdB,f_ML_OA_MM2,'x', 'MarkerEdgeColor', '#A2142F', 'MarkerFaceColor', 'none', 'MarkerSize', 8,'LineWidth', 1.2)
hold on

semilogy(NsePwrVecdB,CRLB,'-k','linewidth',0.5)
legend('WSSF','ML (exhaustive search)','ML-QN','ML-MM (proposed)','CRLB','Fontsize',12)
xlabel('SNR [dB]','Fontsize',12)
ylabel('RMSE [m]','Fontsize',12)

figure(2)
semilogy(NsePwrVecdB(1:2),f_WSSF(1:2), 's', 'MarkerEdgeColor', '#77AC30', 'MarkerFaceColor', 'none', 'MarkerSize', 8, 'LineWidth', 1.2)
hold on
semilogy(NsePwrVecdB(1:2),f_ML_OA(1:2), 'o', 'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor', 'none', 'MarkerSize', 8, 'LineWidth', 1.2)
hold on
semilogy(NsePwrVecdB(1:2),f_ML_OA_MM2(1:2),'x', 'MarkerEdgeColor', '#A2142F', 'MarkerFaceColor', 'none', 'MarkerSize', 8,'LineWidth', 1.2)


figure(3)
yyaxis left
plot(NsePwrVecdB,time_WSSF1*1000,  '-.s', 'Color', '#77AC30', 'MarkerFaceColor', 'none', 'MarkerSize', 8, 'LineWidth', 1.2)
hold on
plot(NsePwrVecdB,time_ML_OA1*1000, '-.o', 'Color', [0 0.4470 0.7410], 'MarkerFaceColor', 'none', 'MarkerSize', 8, 'LineWidth', 1.2)
hold on
plot(NsePwrVecdB,time_GN_OA122*1000,'-.*', 'Color', '#0000FF', 'MarkerSize', 8, 'LineWidth', 1.2)
hold on
plot(NsePwrVecdB,time_ML_OA_MM122*1000,'-.x', 'Color', '#A2142F', 'MarkerSize', 8,'LineWidth', 1.2)
hold on
ylabel('Average Runtime [ms]','Fontsize',16)
yyaxis right
plot(NsePwrVecdB,iteration_ML_GN,'->', 'Color', [0.28 0.49 0.17], 'MarkerSize', 8,'LineWidth', 1.2)
hold on
plot(NsePwrVecdB,iteration_OA_MM,'-^', 'Color', [0.17 0.27 0.86], 'MarkerSize', 8,'LineWidth', 1.2)
ylabel('Number of Iterations','Fontsize',12)
legend('WSSF (runtime)','ML (runtime)','ML-QN (runtime)','ML-MM (runtime)','ML-QN (number of iterations)','ML-MM (number of iterations)','Fontsize',12)
xlabel('SNR [dB]','Fontsize',12)