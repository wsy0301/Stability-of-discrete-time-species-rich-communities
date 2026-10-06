
% This file can be used to generate Fig. 2ab

%% 
clc
clear all


%%
figure
set(gcf,'unit','centimeters','position',[15 5 15 6]) %长11 宽7


%%
subplot(1,2,1)
set(gca,'unit', 'centimeters', 'position', [1.5,1,5,3.5]);
set(gca,'FontSize',8) %坐标轴字体的大小

% A=[-1,0,0.5;
%     0,-1,-0.5;
%     0.5,0.5,-1];

% A=[-1,0.6,0.8;
%     -0.6,-1,0.6;
%     0.6,0.3,-1];

A=[-1,0.5,0.8;
    -0.5,-1,0.6;
    0.6,0.3,-1];

X=ones(3,1);
r=-A*X;

%验证离散中是否成立
diag(X)*exp(r+A*X)

%初始丰度
x_1=1;
x_2=2;
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
legend([h1 h2 h3],{ 'Species 1','Species 2','Species 3'},'fontsize',8)
legend boxoff
box on;
%xlabel('Time','fontsize',9);
set(gca,'YLim',[0 2.3])
%ylabel('Abundance','FontSize',9); 
set(gca,'xminortick','on')  
set(gca,'XLim',[0 17])

%恢复时间
Time(1,1)
%%
subplot(1,2,2)
set(gca,'unit', 'centimeters', 'position', [7,1,5,3.5]);
set(gca,'FontSize',8) %坐标轴字体的大小

% A=[-1,-0.5,0;
%     0.5,-1,-0.5;
%     0,0.5,-1];
A=[-1,-0.6,0.3;
    0.6,-1,0.8;
    0.3,0.6,-1];

X=ones(3,1);
r=-A*X;

%验证离散中是否成立
diag(X)*exp(r+A*X)

%初始丰度
x_1=1;
x_2=2;
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
plot(t,X(1,:),'-','color',[0 0.4470 0.7410],'LineWidth',1);
plot(t,X(2,:),'-','color',[0.8500 0.3250 0.0980],'LineWidth',1);
plot(t,X(3,:),'-','color',[0.9290 0.6940 0.1250],'LineWidth',1);
%
%legend([h_e h_c h_m],{ '+/-','-/-','+/+'},'fontsize',8)
%legend boxoff
box on;
%xlabel('Time','fontsize',9);
set(gca,'YLim',[0 2.3])
%ylabel('Abundance','FontSize',9); 
%set(gca,'xminortick','on')  
set(gca,'XLim',[0 17])

%恢复时间
Time(1,1)

%修改坐标轴刻度距离坐标轴的距离
%tight_subplot(1,1,[.01 .03],[.1 .01],[.01 .01])
