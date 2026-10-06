% This file can be used to generate Fig. 3a


clear all
clc

figure
set(gcf,'unit','centimeters','position',[15 5 12 9]) %长11 宽7
subplot(1,1,1)
set(gca,'unit', 'centimeters', 'position', [2 1 4.5 3]);
set(gca,'FontSize',8) %坐标轴字体的大小


%% 画不稳定的灰色区域
point1=[0 3 3 0];
point2=[-0.4 -0.4 0 0];
Color=[0.75 0.75 0.75];
a=patch(point1,point2,Color) 

alpha(0.6)
set(a,'edgecolor','none')

text(0.25,-0.21,'Unstable','Color','black','FontSize',8)
hold on 
%% 系统参数 随着自调节d的变化 稳定性指标的变化 稳定性指标为log10(max(svd(J_abs)))
syms d
sigma=0.05;

S=100;
C=0.1;

%S=100;
%C=0.2;

%S=100;
%C=0.3;

%S=200;
%C=0.1;

%S=300;  
%C=0.1;

%% 随机
x_lm=-(d-1)-sigma*sqrt(S*C); %最左边特征值的实部
y_lm=0;   %最左边特征值的虚部
x_rm=-(d-1)+sigma*sqrt(S*C); %最右边特征值的实部
y_rm=0; %最右边特征值的虚部
x_um=-(d-1); %最上边特征值的实部
y_um=sigma*sqrt(S*C); %最上边特征值的虚部

abs_lm=sqrt(x_lm^2+y_lm^2);
abs_rm=sqrt(x_rm^2+y_rm^2);
abs_um=sqrt(x_um^2+y_um^2);


i=1;
for val_d=0:0.01:3
val_zlm=double(subs(abs_lm,d,val_d)); %最左边特征值对应的绝对值
val_zrm=double(subs(abs_rm,d,val_d)); %最右边特征值对应的绝对值
val_zum=double(subs(abs_um,d,val_d)); %最上边特征值对应的绝对值
J_abs=[val_zlm,val_zrm,val_zum]; %为了用max 构造矩阵 详情见max用法 
alpha_r{1,i}=-log(max(J_abs)); %求三个特征值对应的绝对值的最大——alpha
i=i+1;
end
val_d=0:0.01:3; 
h_r=plot(val_d,cell2mat(alpha_r),'linestyle','-',"Color",'[0 0.4470 0.7410]','LineWidth',1)
hold on 

%% 生成50个随机矩阵

n=50;
J=cell(1,n);%n为50 生成50个随即交互的群落矩阵

for d=0.1:0.1:3;

for t=1:n;
for i=1:S;
    for j=1:S;
        if i==j;
            M{i,j}=-(d-1);
        else 
             M{i,j}=rand;
            if M{i,j}<=C;
                 M{i,j}=normrnd(0,sigma,1,1);
            else  M{i,j}=0;
            end
        end
    end
end
J{1,t}=M;
end
val_z=cell(1,S);
alpha1=cell(1,n);
for t=1:n
    lam=eig(cell2mat(J{1,t}));%求群落矩阵M的特征值 
for k=1:S
    val_z{1,k}=abs(lam(k)); %求M矩阵的特征值的模
end
alpha1{1,t}=-log(max(cell2mat(val_z)));%每一个群落矩阵M中S个特征值对应的特征根的最大实部——取复数——alpha
end
alpha_mean=mean(cell2mat(alpha1));

h_rr=plot(d,alpha_mean,'d','Markersize',4,'MarkerEdgeColor',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410])
end
axis ([0,2.2,-0.3,2.3])

xlabel('$$s$$','Interpreter','latex','fontsize',9);
ylabel('Stability','FontSize',9); 
box on