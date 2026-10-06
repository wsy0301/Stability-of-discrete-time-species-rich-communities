
% This file can be used to generate Fig. 3b


%% 画示意图
clear all
clc


%% 参数 （随机，真实的+/+ +/-）
sigma=0.05;
S=100;  
C=0.3;
n=5;

%% 画随机相互作用 d=0.5
figure 
set(gcf, 'unit', 'centimeters', 'position', [10 5 23 15]); %设置 Figure 在整个电脑屏幕中的位置及大小
%set(gca,'unit', 'centimeters', 'position', [10, 10, 6, 6]); %设置 Axis(坐标轴）在Figture中的左边界、下边界、宽度、高度
%set(gca,'FontSize',8); %设置坐标轴刻度的字体的大小
%set(gcf,'Position',[100 100 700 150]); %设置 Figure 在整个电脑屏幕中的 位置及大小
%set(gca,'Position',[.13 .17 .750 .72]);%设置 Axis(坐标轴）在Figture中的左边界、下边界、宽度、高度
%figure_FontSize=8;


subplot(3,3,1)
set(gca,'unit', 'centimeters', 'position', [2 11 4 4]);

%画不稳定的区域
point1=[-2 2 2 -2];
point2=[-2 -2 2 2];
Color=[0.75 0.75 0.75];
a=patch(point1,point2,Color) ;
alpha(0.5);
set(a,'edgecolor','none');
hold on ;

%画单位圆
% t=0:0.01:2*pi;
% x=cos(t);
% y=sin(t);
% h_circle=plot(x,y,'linestyle','-',"Color",'black','LineWidth',1.1)

axis equal;

%画圆内表示稳定性的渐变图

r=(0:0.01:1)';
theta=pi*(-1:0.01:1);
x_circle=r*cos(theta);
y_circle=r*sin(theta);

circle=x_circle.^2+y_circle.^2;
p_circle=pcolor(x_circle,y_circle,4.5*(1-circle));
%circle=sqrt(x_circle.^2+y_circle.^2);
%p_circle=pcolor(x_circle,y_circle,-log(circle));
% colormap autumn

% 自定义 colormap 
start_color = [253,174,107]/255; % 开始颜色
end_color = [255,245,235]/255;    % 结束颜色
n_colors = 256;                      % 颜色数量
custom_colormap = [linspace(start_color(1), end_color(1), n_colors)', ... % 红色分量
                   linspace(start_color(2), end_color(2), n_colors)', ... % 绿色分量
                   linspace(start_color(3), end_color(3), n_colors)'];    % 蓝色分量

%应用自定义 colormap
colormap(custom_colormap);



set(p_circle,'LineStyle','none');%去除网格线
%colorbar 
axis([-1.5,1.5,-1.1,1.1]);

xticks([-1 0 1]);
yticks([-1 0 1]);

box on;
xlabel('Real','fontsize',9);
ylabel('Imaginary','FontSize',9);
title('$$s=0.5$$','Interpreter','latex','fontsize',10)



%d=0.5
d=0.5;
J=cell(n,1);
M=cell(S);
for t=1:n
for i=1:S
    for j=1:S
        if i==j
            M{i,j}=-(d-1);
        else 
             M{i,j}=rand;
            if M{i,j}<=C
                 M{i,j}=normrnd(0,sigma,1,1);
            else  M{i,j}=0;
            end
        end
    end
end
J{t,1}=M;
end
J51=cell2mat(J{5,1});
J51_lam=eig(J51);
h51=plot(real(J51_lam),imag(J51_lam),'.',"MarkerEdgeColor",[31,120,180]/255)

%画随机相互作用的圆
t=linspace(0,2*pi,1000);
x_r=-(d-1)+sigma*sqrt(S*C)*cos(t);
y_r=sigma*sqrt(S*C)*sin(t);
h_r=plot(x_r,y_r,'linestyle','-',"Color",'black','LineWidth',1)

%% 画随机相互作用 d=1

subplot(3,3,2)
set(gca,'unit', 'centimeters', 'position', [7 11 4 4]);


%画不稳定的区域
point1=[-2 2 2 -2];
point2=[-2 -2 2 2];
Color=[0.75 0.75 0.75];
a=patch(point1,point2,Color) ;
alpha(0.5);
set(a,'edgecolor','none');
hold on ;

%画单位圆
% t=0:0.01:2*pi;
% x=cos(t);
% y=sin(t);
% h_circle=plot(x,y,'linestyle','-',"Color",'black','LineWidth',1.1)

axis equal;

%画圆内表示稳定性的渐变图
r=(0:0.01:1)';
theta=pi*(-1:0.01:1);
x_circle=r*cos(theta);
y_circle=r*sin(theta);

circle=x_circle.^2+y_circle.^2;
p_circle=pcolor(x_circle,y_circle,4.5*(1-circle));
% 自定义 colormap 
start_color = [253,174,107]/255; % 开始颜色
end_color = [255,245,235]/255;    % 结束颜色
n_colors = 256;                      % 颜色数量
custom_colormap = [linspace(start_color(1), end_color(1), n_colors)', ... % 红色分量
                   linspace(start_color(2), end_color(2), n_colors)', ... % 绿色分量
                   linspace(start_color(3), end_color(3), n_colors)'];    % 蓝色分量

%应用自定义 colormap
colormap(custom_colormap);
%colormap summer
set(p_circle,'LineStyle','none');%去除网格线
 
axis([-1.5,1.5,-1.1,1.1]);

xticks([-1 0 1]);
yticks([-1 0 1]);

box on;
title('$$s=1$$','Interpreter','latex','fontsize',10);

%d=1
d=1;
J=cell(n,1);
M=cell(S);
for t=1:n
for i=1:S
    for j=1:S
        if i==j
            M{i,j}=-(d-1);
        else 
             M{i,j}=rand;
            if M{i,j}<=C
                 M{i,j}=normrnd(0,sigma,1,1);
            else  M{i,j}=0;
            end
        end
    end
end
J{t,1}=M;
end

J51=cell2mat(J{5,1});
J51_lam=eig(J51);
h51=plot(real(J51_lam),imag(J51_lam),'.',"MarkerEdgeColor",[31,120,180]/255)%[0.6350 0.0780 0.1840])

%画随机相互作用的圆
t=linspace(0,2*pi,1000);
x_r=-(d-1)+sigma*sqrt(S*C)*cos(t);
y_r=sigma*sqrt(S*C)*sin(t);
h_r=plot(x_r,y_r,'linestyle','-',"Color",'black','LineWidth',1)

%% 画随机相互作用 d=1.5

subplot(3,3,3)
set(gca,'unit', 'centimeters', 'position', [12 11 4 4]);

%画不稳定的区域
point1=[-2 2 2 -2];
point2=[-2 -2 2 2];
Color=[0.75 0.75 0.75];
a=patch(point1,point2,Color) ;
alpha(0.5);
set(a,'edgecolor','none');
hold on ;

%画单位圆
% t=0:0.01:2*pi;
% x=cos(t);
% y=sin(t);
% h_circle=plot(x,y,'linestyle','-',"Color",'black','LineWidth',1.1)

axis equal;

%画圆内表示稳定性的渐变图
r=(0:0.01:1)';
theta=pi*(-1:0.01:1);
x_circle=r*cos(theta);
y_circle=r*sin(theta);

circle=x_circle.^2+y_circle.^2;
p_circle=pcolor(x_circle,y_circle,4.5*(1-circle));

% 自定义 colormap 
start_color = [253,174,107]/255; % 开始颜色
end_color = [255,245,235]/255;    % 结束颜色
n_colors = 256;                      % 颜色数量
custom_colormap = [linspace(start_color(1), end_color(1), n_colors)', ... % 红色分量
                   linspace(start_color(2), end_color(2), n_colors)', ... % 绿色分量
                   linspace(start_color(3), end_color(3), n_colors)'];    % 蓝色分量

%应用自定义 colormap
colormap(custom_colormap);


set(p_circle,'LineStyle','none');%去除网格线

axis([-1.5,1.5,-1.1,1.1]);

xticks([-1 0 1]);
yticks([-1 0 1]);

box on;
title('$$s=1.5$$','Interpreter','latex','fontsize',10);

%d=1.5
d=1.5;
J=cell(n,1);
M=cell(S);
for t=1:n
for i=1:S
    for j=1:S
        if i==j
            M{i,j}=-(d-1);
        else 
             M{i,j}=rand;
            if M{i,j}<=C
                 M{i,j}=normrnd(0,sigma,1,1);
            else  M{i,j}=0;
            end
        end
    end
end
J{t,1}=M;
end

J51=cell2mat(J{5,1});
J51_lam=eig(J51);
h51=plot(real(J51_lam),imag(J51_lam),'.',"MarkerEdgeColor",[31,120,180]/255)%[0.6350 0.0780 0.1840])


%画随机相互作用的圆
t=linspace(0,2*pi,1000);
x_r=-(d-1)+sigma*sqrt(S*C)*cos(t);
y_r=sigma*sqrt(S*C)*sin(t);
h_r=plot(x_r,y_r,'linestyle','-',"Color",'black','LineWidth',1)

%% 
subplot(3,3,9)
set(gca,'unit', 'centimeters', 'position', [10 3 3 3]);
%画不稳定的区域
point1=[-2 2 2 -2];
point2=[-2 -2 2 2];
Color=[0.75 0.75 0.75];
a=patch(point1,point2,Color) ;
alpha(0.5);
set(a,'edgecolor','none');
hold on ;

%画单位圆
% t=0:0.01:2*pi;
% x=cos(t);
% y=sin(t);
% h_circle=plot(x,y,'linestyle','-',"Color",'black','LineWidth',1.1)
axis equal;

%画圆内表示稳定性的渐变图
r=(0:0.01:1)';
theta=pi*(-1:0.01:1);
x_circle=r*cos(theta);
y_circle=r*sin(theta);

circle=x_circle.^2+y_circle.^2;
p_circle=pcolor(x_circle,y_circle,4.5*(1-circle));


% 自定义 colormap 
start_color = [253,174,107]/255; % 开始颜色
end_color = [255,245,235]/255;    % 结束颜色
n_colors = 256;                      % 颜色数量
custom_colormap = [linspace(start_color(1), end_color(1), n_colors)', ... % 红色分量
                   linspace(start_color(2), end_color(2), n_colors)', ... % 绿色分量
                   linspace(start_color(3), end_color(3), n_colors)'];    % 蓝色分量

%应用自定义 colormap
colormap(custom_colormap);
%
c=colorbar;
%
set(c,'position', [.8 .2 .01 .2]) %左右位置 上下位置 宽度 高度

set(p_circle,'LineStyle','none');%去除网格线

axis ([-1.1,1.1,-1.1,1.1]);
box on;

 