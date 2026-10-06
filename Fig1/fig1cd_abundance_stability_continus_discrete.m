
% This file can be used to generate Figs. 1c-d

clear all
clc

figure;
set(gcf, 'unit', 'centimeters', 'position', [10 5 20 10]);

%% 不稳定的 连续的 丰度响应图
% subplot(2,4,1)
% set(gca,'unit', 'centimeters', 'position', [2,6,3.5,2.5]); %设置 Axis(坐标轴）在Figture中的左边界、下边界、宽度、高度
% set(gca,'FontSize',8)
% 
% %不稳定的
% %A_value=[-1,0.95,0;1,-1,0.95;0,0.95,-1];
% %A_value=[-1,0.6,0.3;0.6,-1,0.8;0.3,0.8,-1];
% 
% % A_value=[-1,-1,0;
% %          -1,-1,1;
% %          0,1,-1];
% 
% A_value=[-1,-1,0;
%          -1,-1,1;
%          0,1,-1];
% 
% X_e=ones(3,1);
% r_value=-A_value*X_e
% 
% %初始丰度
% x_1=2;
% x_2=1;
% x_3=1;
% x_0=[x_1;x_2;x_3];
% 
% 
% 
% %求解
% t=[1:0.01:200];
% %[T,X]=ode45('odefun',t,x_0)
% options=odeset('RelTol',1e-10,'AbsTol',1e-15) 
% [T,X]=ode45('odefun',t,x_0,options);
% 
% 
% %画图
% hold on
% h1=plot(T,X(:,1),'-','color',[0 0.4470 0.7410],'LineWidth',1);
% h2=plot(T,X(:,2),'-','color',[0.8500 0.3250 0.0980],'LineWidth',1);
% h3=plot(T,X(:,3),'-','color',[0.9290 0.6940 0.1250],'LineWidth',1);
%  
%  box on;
%  set(gca,'YLim',[0 3.3]);
% % set(gca,'xminortick','on')  
%  set(gca,'XLim',[0 15.3]);
%  
%  ylabel('Abundance','FontSize',9); 
% xlabel('Generation','FontSize',9); 
% title('Response','FontSize',10);
% 
% %% 不稳定的 连续的 特征值分布
% %群落矩阵M、连接矩阵A、物种丰度X在平衡点的关系为：M=diag{X}A
% %内禀增长率可根据平衡点定义f=0 逆求出来
% %因此，我们考虑平衡点X=1 直接给定矩阵M 
% %分别画出稳定的、不稳定的M的特征值分布
% 
% subplot(2,4,2)
% set(gca,'unit', 'centimeters', 'position', [6.5,6,3.5,2.5]); %设置 Axis(坐标轴）在Figture中的左边界、下边界、宽度、高度
% set(gca,'FontSize',8)
% 
% hold on
% xline (0,'k','LineWidth',0.5);
% 
% %画不稳定的区域
% point1=[-3 3 3 -3];
% point2=[-3 -3 3 3];
% Color=[0.75 0.75 0.75];
% a=patch(point1,point2,Color) ;
% alpha(0.4);
% set(a,'edgecolor','none');
% 
% 
% %画稳定的左半平面
% point1=[-3 0 0 -3];
% point2=[-3 -3 3 3];
% %Color=[204 222 216]/255;
% Color=[200 255 239]/255;
% a=patch(point1,point2,Color) ;
% %alpha(0.6);
% set(a,'edgecolor','none');
% hold on;
% 
% 
% 
% %不稳定的
% %M_unstable=[-1,0.95,0;1,-1,0.95;0,0.95,-1];
% %M_unstable=[-1,0.6,0.3;0.6,-1,0.8;0.3,0.8,-1];
% M_unstable=A_value;
% eig_unstable=eig(M_unstable);
% h_unstable=plot(real(eig_unstable),imag(eig_unstable),'.',"MarkerEdgeColor",'k',"MarkerSize",10)
% 
% axis ([-2.6,1.48,-1.1,1.1]);
% 
% 
% box on
% 
% xlabel('Real','fontsize',9);
% ylabel('Imaginary','fontsize',9); 
% title('Stability criterion','FontSize',10);
% 
% %% 离散的 不稳定的 物种丰度的响应
% subplot(2,4,3)
% set(gca,'unit', 'centimeters', 'position', [11,6,3.5,2.5]); %设置 Axis(坐标轴）在Figture中的左边界、下边界、宽度、高度
% set(gca,'FontSize',8)
% 
% 
% X=ones(3,1);
% %A=[-1,0.6,0.3;0.6,-1,0.8;0.3,0.8,-1];
% A=A_value;
% r=-A*X;
% 
% 
% %初始丰度
% x_1=2;
% x_2=1;
% x_3=1;
% X_0=[x_1;x_2;x_3];
% 
% Y{1,1}=X_0;
% t=1;
% for i=2:1:25;
% Y{1,i}=diag(Y{1,i-1})*exp(r+A*Y{1,i-1});
% %if abs(cell2mat(Y(1,i))-ones(3,1))<0.001*ones(3,1) & abs(cell2mat(Y(1,i-1))-ones(3,1))<0.001*ones(3,1)
% %    Time(t,1)=i-1;
% %    t=t+1;
% %end
% end
%  
% X=cell2mat(Y);
% hold on 
% t=1:1:25;
% h1=plot(t,X(1,:),'-','color',[0 0.4470 0.7410],'LineWidth',1);
% h2=plot(t,X(2,:),'-','color',[0.8500 0.3250 0.0980],'LineWidth',1);
% h3=plot(t,X(3,:),'-','color',[0.9290 0.6940 0.1250],'LineWidth',1);
% 
% %
% box on;
% 
% set(gca,'YLim',[0 3.3]);
% %set(gca,'xminortick','on')  
% set(gca,'XLim',[0 8.3]);
% 
% ylabel('Abundance','FontSize',9); 
% xlabel('Generation','FontSize',9); 
% title('Response','FontSize',10);
% 
% %% 离散的 不稳定的 特征值分布
% subplot(2,4,4)
% set(gca,'unit', 'centimeters', 'position', [16,6,3.5,2.5]); %设置 Axis(坐标轴）在Figture中的左边界、下边界、宽度、高度
% set(gca,'FontSize',8)
% 
% 
% %画不稳定的区域
% point1=[-3 3 3 -3];
% point2=[-3 -3 3 3];
% Color=[0.75 0.75 0.75];
% a=patch(point1,point2,Color) ;
% alpha(0.5);
% set(a,'edgecolor','none');
% hold on;
% 
% %画单位圆
% t=0:0.01:2*pi;
% x=cos(t);
% y=sin(t);
% h_circle=plot(x,y,'linestyle','-',"Color",'black','LineWidth',1)
% axis equal;
% 
% %画圆内表示稳定性的渐变图
% r=(0:0.01:1)';
% theta=pi*(-1:0.01:1);
% x_circle=r*cos(theta);
% y_circle=r*sin(theta);
% 
% circle=x_circle.^2+y_circle.^2;
% circle_color=ones(size(circle));
% p_circle=pcolor(x_circle,y_circle,circle_color);
% 
% map=[200 255 239]/255;
% colormap(map)
% 
% set(p_circle,'LineStyle','none');%去除网格线
% 
% 
% 
% % M=I+diag{X}A
% %A_unstable=[-1,0.6,0.3;0.6,-1,0.8;0.3,0.8,-1];
% A_unstable=A_value;
% M_unstable=eye(3)+A_unstable;
% 
% eig_unstable=eig(M_unstable);
% h_unstable=plot(real(eig_unstable),imag(eig_unstable),'.',"MarkerEdgeColor",'k',"MarkerSize",10)
% 
% axis ([-1.54,1.54,-1.1,1.1]);
% 
% box on
% xlabel('Real','fontsize',9);
% ylabel('Imaginary','fontsize',9); 
% title('Stability criterion','FontSize',10);

%% 稳定的 连续的 物种丰度的响应图
subplot(2,4,5)
set(gca,'unit', 'centimeters', 'position', [2,1,3.5,2.5]); %设置 Axis(坐标轴）在Figture中的左边界、下边界、宽度、高度
set(gca,'FontSize',8)

%稳定的
%A_value=[-1,-1,0;1,-1,-0.5;0,0.5,-1];
A_value=[-1,-0.3,0.6;
        -0.3,-1,0.7;
        0.6,0.7,-1];
X_e=ones(3,1);
r_value=-A_value*X_e;

%初始丰度
x_1=2;
x_2=1;
x_3=1;
x_0=[x_1;x_2;x_3];

%
hold on
t=[1:0.01:20];
%[T,X]=ode45('odefun',t,x_0);
options=odeset('RelTol',1e-20,'AbsTol',1e-25) 
[T,X]=ode45('odefun',t,x_0,options);

h1=plot(t,X(:,1),'-','color',[0 0.4470 0.7410],'LineWidth',1);
h2=plot(t,X(:,2),'-','color',[0.8500 0.3250 0.0980],'LineWidth',1);
h3=plot(t,X(:,3),'-','color',[0.9290 0.6940 0.1250],'LineWidth',1);

%
legend([h1 h2 h3],{ 'Species 1','Species 2','Species 3'},'fontsize',8);
legend boxoff;
box on;
set(gca,'YLim',[0 3.3]);
%set(gca,'xminortick','on')
set(gca,'XLim',[0 15.3]);
ylabel('Abundance','FontSize',9); 
xlabel('Generation','FontSize',9); 
title('Response','FontSize',10);


%% 连续的 稳定的 特征值分布
subplot(2,4,6)
set(gca,'unit', 'centimeters', 'position', [6.5,1,3.5,2.5]); %设置 Axis(坐标轴）在Figture中的左边界、下边界、宽度、高度
set(gca,'FontSize',8)

hold on
xline (0,'k','LineWidth',0.5);

%画不稳定的区域
point1=[-3 3 3 -3];
point2=[-3 -3 3 3];
Color=[0.75 0.75 0.75];
a=patch(point1,point2,Color) ;
alpha(0.5);
set(a,'edgecolor','none');


%画稳定的左半平面
%画不稳定的区域
point1=[-3 0 0 -3];
point2=[-3 -3 3 3];
Color=[200 255 239]/255;
a=patch(point1,point2,Color) ;
%alpha(0.6);
set(a,'edgecolor','none');
hold on;


%稳定的
%M_stable=[-1,-1,0;1,-1,-0.5;0,0.5,-1];
%M_stable=[-1,-0.3,0.6;-0.3,-1,0.7;0.6,0.7,-1];
M_stable=A_value;
eig_stable=eig(M_stable);
h_stable=plot(real(eig_stable),imag(eig_stable),'.',"MarkerEdgeColor",'k',"MarkerSize",10)

axis ([-2.6,1.48,-1.1,1.1]);

%ylim([-1.1 1.1])


box on
xlabel('Real','fontsize',9);
ylabel('Imaginary','fontsize',9); 
title('Stability criterion','FontSize',10);



%% 离散的 稳定的 物种丰度的响应
subplot(2,4,7)
set(gca,'unit', 'centimeters', 'position', [11,1,3.5,2.5]); %设置 Axis(坐标轴）在Figture中的左边界、下边界、宽度、高度
set(gca,'FontSize',8)

X_e=ones(3,1);
%A=[-1,0.3,0;0.3,-1,0.3;0,0.3,-1];
A=A_value;
r=-A*X_e;
% 验证是否X_e=f(X_e)
diag(X_e)*exp(r+A*X_e);

%初始丰度
x_1=2;
x_2=1;
x_3=1;
X_0=[x_1;x_2;x_3];

Y{1,1}=X_0;
t=1;
for i=2:1:25;
Y{1,i}=diag(Y{1,i-1})*exp(r+A*Y{1,i-1});
if abs(cell2mat(Y(1,i))-ones(3,1))<0.001*ones(3,1) & abs(cell2mat(Y(1,i-1))-ones(3,1))<0.001*ones(3,1)
    Time(t,1)=i-1;
    t=t+1;
end
end
 
X=cell2mat(Y);
hold on 
t=1:1:25;
h1=plot(t,X(1,:),'-','color',[0 0.4470 0.7410],'LineWidth',1);
h2=plot(t,X(2,:),'-','color',[0.8500 0.3250 0.0980],'LineWidth',1);
h3=plot(t,X(3,:),'-','color',[0.9290 0.6940 0.1250],'LineWidth',1);

%
box on;
set(gca,'YLim',[0 3.3]);
%set(gca,'xminortick','on')  
set(gca,'XLim',[0 8.3]);

ylabel('Abundance','FontSize',9); 
xlabel('Generation','FontSize',9); 
title('Response','FontSize',10);

%% 离散的 稳定的  特征值分布
subplot(2,4,8)
set(gca,'unit', 'centimeters', 'position', [16,1,3.5,2.5]); %设置 Axis(坐标轴）在Figture中的左边界、下边界、宽度、高度
set(gca,'FontSize',8)

%画不稳定的区域
point1=[-3 3 3 -3];
point2=[-3 -3 3 3];
Color=[0.75 0.75 0.75];
a=patch(point1,point2,Color) ;
alpha(0.5);
set(a,'edgecolor','none');
hold on;

%画单位圆
t=0:0.01:2*pi;
x=cos(t);
y=sin(t);
h_circle=plot(x,y,'linestyle','-',"Color",'black','LineWidth',1)
axis equal;

%画圆内表示稳定性的渐变图
r=(0:0.01:1)';
theta=pi*(-1:0.01:1);
x_circle=r*cos(theta);
y_circle=r*sin(theta);

circle=x_circle.^2+y_circle.^2;
circle_color=ones(size(circle));
p_circle=pcolor(x_circle,y_circle,circle_color);

map=[200 255 239]/255;
colormap(map)
set(p_circle,'LineStyle','none');%去除网格线


% M=I+diag{X}A
%A_stable=[-1,0.3,0;0.3,-1,0.3;0,0.3,-1];
A_stable=A_value;
M_stable=eye(3)+A_stable;

eig_stable=eig(M_stable);
h_stable=plot(real(eig_stable),imag(eig_stable),'.',"MarkerEdgeColor",'k',"MarkerSize",10)

axis ([-1.54,1.54,-1.1,1.1]);

box on
xlabel('Real','fontsize',9);
ylabel('Imaginary','fontsize',9); 
title('Stability criterion','FontSize',10);

