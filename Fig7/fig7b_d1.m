
% This file can be used to generate Fig. 7b
% d=1



clear all
clc


figure;
set(gcf, 'unit', 'centimeters', 'position', [10 5 20 6]); %设置 Figure 在整个电脑屏幕中的位置及大小


S=200;
d=1;
sigma_value=0.1;
C_value=0.1;

n=100;%每个点计算的个数

%% 无结构的捕食-被捕食

subplot(1,3,1)
set(gca,'unit', 'centimeters', 'position', [1.5, 1, 5, 3.5]); %设置 Axis(坐标轴）在Figture中的左边界、下边界、宽度、高度
set(gca,'FontSize',8);
hold on 

%C变化
%C_uns=0.025:0.025:0.5;%无结构的捕食被捕食C的取值范围

C_uns=0.01:0.025:0.4;
sigma=sigma_value;%固定sigma

m=1;
for C=C_uns
    
J_uns=zeros(n,1);
for t=1:n
    M=zeros(S);
for i=1:S
    for j=1:S
        if i==j
            M(i,j)=1-d;
        elseif i<j
                p1=rand;
            if p1<=C
                p2=rand;
                if p2<=0.5
                    M(i,j)=abs(normrnd(0,sigma,1,1));
                    M(j,i)=-abs(normrnd(0,sigma,1,1));
                else 
                    M(i,j)=-abs(normrnd(0,sigma,1,1));
                    M(j,i)=abs(normrnd(0,sigma,1,1));
                end
           end
        end
    end
end
J_uns(t,1)=-log(max(abs(eig(M))));
end
J_1=J_uns;
J_1(J_1>0)=1;
J_1(J_1<=0)=0;
J_unsC(1,m)=sum(J_1,"all")/n;
m=m+1;
end
%画图
B_Cuns=sigma*sqrt(S*C_uns);
h_uns=plot(B_Cuns,J_unsC,'-.x',"Color",'k','LineWidth',1,'Markersize',5)

%临界值
E_absnormal=sqrt(2/pi)*sigma;
syms C_normal
eqn=sigma*sqrt(S*C_normal)*(1+E_absnormal^2/(sigma^2))==1;
C_normal=solve(eqn,C_normal);
hold on 
xline(double(sigma*sqrt(S*vpa(C_normal))),'-','LineWidth',1,'Color',[0.8500 0.3250 0.0980]);


%sigma 变化
%sigma_uns=0.01:0.01:0.2;%无结构的捕食被捕食者的sigma的取值范围

sigma_uns=0.01:0.01:0.15;
C=C_value;%固定C
m=1;
for sigma=sigma_uns
    
J_uns=zeros(n,1);
for t=1:n
    M=zeros(S);
for i=1:S
    for j=1:S
        if i==j
            M(i,j)=1-d;
        elseif i<j
                p1=rand;
            if p1<=C
                p2=rand;
                if p2<=0.5
                    M(i,j)=abs(normrnd(0,sigma,1,1));
                    M(j,i)=-abs(normrnd(0,sigma,1,1));
                else 
                    M(i,j)=-abs(normrnd(0,sigma,1,1));
                    M(j,i)=abs(normrnd(0,sigma,1,1));
                end
           end
        end
    end
end
J_uns(t,1)=-log(max(abs(eig(M))));
end
J_1=J_uns;
J_1(J_1>0)=1;
J_1(J_1<=0)=0;
J_unssigma(1,m)=sum(J_1,"all")/n;
m=m+1;
end

B_sigmauns=sigma_uns*sqrt(S*C);
h_sigmauns=plot(B_sigmauns,J_unssigma,':o',"Color",'k','LineWidth',1,'Markersize',3)

syms sigma_normal
E_absnormal=sqrt(2/pi)*sigma_normal;
eqn=sigma_normal*sqrt(S*C)*(1+E_absnormal^2/(sigma_normal^2))==1;
sigma_normal=solve(eqn,sigma_normal);
xline(double(vpa(sigma_normal)*sqrt(S*C)),'-','LineWidth',1,'Color',[0.8500 0.3250 0.0980]);

%
%xlabel('$$\sigma\sqrt{SC}$$','Interpreter','latex','fontsize',10);
ylabel('$$P(\emph{\emph{Stable}})$$','Interpreter','latex','fontsize',9); 
%axis ([0.04,0.85,-0.1,1.1])
axis ([0,0.93,-0.1,1.1])

%title('Unstructured exploitative (+/-)','FontSize',10)
title('Unstructured food web','FontSize',10)
%text(0.58,1.05,'Unstructured Predator-Prey','Color','black','FontSize',8)
box on

legend([h_uns h_sigmauns],{'$$C$$', '$$\sigma$$'},'Interpreter','latex','fontsize',8)
legend boxoff;%图例去边框

%% 级联捕食被捕食
subplot(1,3,2)
set(gca,'unit', 'centimeters', 'position', [8, 1, 5, 3.5]); %设置 Axis(坐标轴）在Figture中的左边界、下边界、宽度、高度
set(gca,'FontSize',8);

%无结构的捕食-被捕食
hold on
plot(B_Cuns,J_unsC,'-.',"Color",'k','LineWidth',1);
plot(B_sigmauns,J_unssigma,':',"Color",'k','LineWidth',1);
%无结构的捕食-被捕食的临界值
xline(double(vpa(sigma_normal)*sqrt(S*C)),'-','LineWidth',1,'Color',[0.8500 0.3250 0.0980]);

%改变C
%C_cas=[0.0125:0.0125:0.1,0.125:0.025:0.5];%级联的捕食-被捕食的C的取值范围

C_cas=[0.0125:0.0125:0.1,0.125:0.025:0.4];
sigma=sigma_value;%确定sigma
m=1;
for C=C_cas
    
J_cas=zeros(n,1);
for t=1:n
    M=zeros(S);
for i=1:S
    for j=1:S
        if i==j
            M(i,j)=1-d;
        elseif i<j
                p1=rand;
            if p1<=C
                    M(i,j)=-abs(normrnd(0,sigma,1,1));
                    M(j,i)=abs(normrnd(0,sigma,1,1));
            end
        end
    end
end
J_cas(t,1)=-log(max(abs(eig(M))));
end
J_2=J_cas;
J_2(J_2>0)=1;
J_2(J_2<=0)=0;
J_casC(1,m)=sum(J_2,"all")/n;
m=m+1;
end
%画图
B_Ccas=sigma*sqrt(S*C_cas);
h_cas=plot(B_Ccas,J_casC,'-.x',"Color",[0 0.4470 0.7410],'LineWidth',1,'Markersize',5)
hold on

%改变sigma
%sigma_cas=0.01:0.01:0.2;%sigma的取值范围
sigma_cas=0.01:0.01:0.15;
C=C_value;
m=1;
for sigma=sigma_cas
    
J_cas=zeros(n,1);
for t=1:n
    M=zeros(S);
for i=1:S
    for j=1:S
        if i==j
            M(i,j)=1-d;
        elseif i<j
                p1=rand;
            if p1<=C
                    M(i,j)=-abs(normrnd(0,sigma,1,1));
                    M(j,i)=abs(normrnd(0,sigma,1,1));
            end
        end
    end
end
J_cas(t,1)=-log(max(abs(eig(M))));
end
J_2=J_cas;
J_2(J_2>0)=1;
J_2(J_2<=0)=0;
J_cassigma(1,m)=sum(J_2,"all")/n;
m=m+1;
end
%画图
B_sigmacas=sigma_cas*sqrt(S*C);
h_sigmacas=plot(B_sigmacas,J_cassigma,':o',"Color",[0 0.4470 0.7410],'LineWidth',1,'Markersize',3)

%
xlabel('$$\sigma\sqrt{SC}$$','Interpreter','latex','fontsize',9);
%axis ([0.04,0.85,-0.1,1.1])
axis ([0,0.93,-0.1,1.1])
title('Cascade food web','FontSize',10)

%text(0.64,1.05,'Cascade Predator-Prey','Color','black','FontSize',8)
box on

legend([h_cas h_sigmacas],{'$$C$$', '$$\sigma$$'},'Interpreter','latex','fontsize',8)
legend boxoff;%图例去边框
%% 生态位捕食-被捕食
subplot(1,3,3)
set(gca,'unit', 'centimeters', 'position', [14.5, 1, 5, 3.5]); %设置 Axis(坐标轴）在Figture中的左边界、下边界、宽度、高度
set(gca,'FontSize',8);

hold on

%无结构的捕食-被捕食
hold on
plot(B_Cuns,J_unsC,'-.',"Color",'k','LineWidth',1);
plot(B_sigmauns,J_unssigma,':',"Color",'k','LineWidth',1);
%无结构的捕食-被捕食的临界值
xline(double(vpa(sigma_normal)*sqrt(S*C)),'-','LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
%级联的捕食被捕食
plot(B_Ccas,J_casC,'-.',"Color",[0 0.4470 0.7410],'LineWidth',1);
plot(B_sigmacas,J_cassigma,':',"Color",[0 0.4470 0.7410],'LineWidth',1);

%C变化
%C_nic=[0.0125:0.0125:0.1,0.125:0.025:0.5];

C_nic=[0.0125:0.0125:0.1,0.125:0.025:0.4];
sigma=sigma_value;
m=1;
for C=C_nic
    
J_nic=zeros(n,1);
for t=1:n
    
eta=zeros(S,1);    
for i=1:S %生成S个物种的生态位值，生态位半径，生态位中心
    eta(i,1)=rand;
end
eta=sort(eta);

b=zeros(S,1);
r=zeros(S,1);
c=zeros(S,1);
for i=1:S
    b(i,1)=betarnd(1,1/C-1);
    r(i,1)=eta(i,1)*b(i,1);
    c(i,1)=r(i,1)/2+(min(1-r(i,1)/2,eta(i,1))-r(i,1)/2)*rand;
end %生成生态位结束

K=zeros(S);
for j=1:S %判断物种与物种之间的生态位 确定邻接矩阵K
    for i=1:S
        if (eta(j,1)>=(c(i,1)-r(i,1)/2) && eta(j,1)<=(c(i,1)+r(i,1)/2))
            K(i,j)=1;
        end
    end
end %生成K结束

M=zeros(S);
for i=1:S %根据邻接矩阵K 生成群落矩阵M
    for j=1:S
        if  i~=j 
            if (K(i,j)==1 && K(j,i)==1)
            p=rand;
                if p<=0.5
                    M(i,j)=abs(normrnd(0,sigma,1,1));
                    M(j,i)=-abs(normrnd(0,sigma,1,1));
                else 
                    M(i,j)=-abs(normrnd(0,sigma,1,1));
                    M(j,i)=abs(normrnd(0,sigma,1,1));
                end
        elseif  (K(i,j)==1 && K(j,i)==0)
            M(i,j)=abs(normrnd(0,sigma,1,1));
            M(j,i)=-abs(normrnd(0,sigma,1,1));
            end
        end
    end
end %生成M结束
M(logical(eye(size(M))))=(1-d)*ones(S,1);
J_nic(t,1)=-log(max(abs(eig(M))));
end
J_3=J_nic;
J_3(J_3>0)=1;
J_3(J_3<=0)=0;
J_nicC(1,m)=sum(J_3,"all")/n;
m=m+1;
end
%画图
B_Cnic=sigma*sqrt(S*C_nic);
h_nic=plot(B_Cnic,J_nicC,'-.x',"Color",[0.9290 0.6940 0.1250],'LineWidth',1,'Markersize',5)

%sigma 变化
%sigma_nic=0.01:0.01:0.2;

sigma_nic=0.01:0.01:0.15;
C=C_value;
m=1;
for sigma=sigma_nic
    
J_nic=zeros(n,1);

for t=1:n
    
eta=zeros(S,1);    
for i=1:S %生成S个物种的生态位值，生态位半径，生态位中心
    eta(i,1)=rand;
end
eta=sort(eta);

b=zeros(S,1);
r=zeros(S,1);
c=zeros(S,1);
for i=1:S
    b(i,1)=betarnd(1,1/C-1);
    r(i,1)=eta(i,1)*b(i,1);
    c(i,1)=r(i,1)/2+(min(1-r(i,1)/2,eta(i,1))-r(i,1)/2)*rand;
end %生成生态位结束

K=zeros(S);
for j=1:S %判断物种与物种之间的生态位 确定邻接矩阵K
    for i=1:S
        if (eta(j,1)>=(c(i,1)-r(i,1)/2)&& eta(j,1)<=(c(i,1)+r(i,1)/2))
            K(i,j)=1;
        end
    end
end %生成K结束

M=zeros(S);
for i=1:S %根据邻接矩阵K 生成群落矩阵M
    for j=1:S
        if  i~=j 
            if (K(i,j)==1 && K(j,i)==1)
            p=rand;
                if p<=0.5
                    M(i,j)=abs(normrnd(0,sigma,1,1));
                    M(j,i)=-abs(normrnd(0,sigma,1,1));
                else 
                    M(i,j)=-abs(normrnd(0,sigma,1,1));
                    M(j,i)=abs(normrnd(0,sigma,1,1));
                end
        elseif  (K(i,j)==1 && K(j,i)==0)
            M(i,j)=abs(normrnd(0,sigma,1,1));
            M(j,i)=-abs(normrnd(0,sigma,1,1));
            end
        end
    end
end %生成M结束
M(logical(eye(size(M))))=(1-d)*ones(S,1);
J_nic(t,1)=-log(max(abs(eig(M))));
end
J_3=J_nic;
J_3(J_3>0)=1;
J_3(J_3<=0)=0;
J_nicsigma(1,m)=sum(J_3,"all")/n;
m=m+1;
end
%画图
B_sigmanic=sigma_nic*sqrt(S*C);
h_sigmanic=plot(B_sigmanic,J_nicsigma,':o',"Color",[0.9290 0.6940 0.1250],'LineWidth',1,'Markersize',3)

%axis ([0.04,0.85,-0.1,1.1])
axis ([0,0.93,-0.1,1.1])
title('Niche food web','FontSize',10)
%text(0.69,1.05,'Niche Predator-Prey','Color','black','FontSize',8)
box on

legend([h_nic h_sigmanic],{'$$C$$', '$$\sigma$$'},'Interpreter','latex','fontsize',8)
legend boxoff;%图例去边框