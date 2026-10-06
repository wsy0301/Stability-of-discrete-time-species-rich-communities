

% This file can be used to generate Fig. 2ab

%% 
clc
clear all


%%
figure
set(gcf,'unit','centimeters','position',[15 5 15 6]) %长11 宽7


%%
subplot(1,2,1)
set(gca,'unit', 'centimeters', 'position', [1.5, 1, 5, 3.5]);
set(gca,'FontSize',8) %坐标轴字体的大小

%画不稳定的区域
point1=[-3 3 3 -3];
point2=[-3 -3 3 3];
Color=[0.75 0.75 0.75];
a=patch(point1,point2,Color) ;
alpha(0.6);
set(a,'edgecolor','none');
hold on;

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
%colorbar 


%画不同半径的稳定性
% r=0.8;
r1=0.8;
theta=pi*(-1:0.01:1);
x_1=r1*cos(theta);
y_1=r1*sin(theta);
h_1=plot(x_1,y_1,'linestyle','--',"Color",[82 82 82]/255,'LineWidth',1)
stability_1=-log(r1)

%r=0.6;
r2=0.6;
theta=pi*(-1:0.01:1);
x_2=r2*cos(theta);
y_2=r2*sin(theta);
h_2=plot(x_2,y_2,'linestyle',':',"Color",[82 82 82]/255,'LineWidth',1)
stability_2=-log(r2)

%r=0.4;
r3=0.4;
theta=pi*(-1:0.01:1);
x_3=r3*cos(theta);
y_3=r3*sin(theta);
h_3=plot(x_3,y_3,'linestyle','-',"Color",[82 82 82]/255,'LineWidth',1)
stability_3=-log(r3)




%画群落矩阵的特征值分布
A=[-1,0.5,0.8;
    -0.5,-1,0.6;
    0.6,0.3,-1];

M=eye(3)+A;

M_eig=eig(M)

%hM_1=plot(real(M_eig),imag(M_eig),'o',"MarkerEdgeColor",'black',"MarkerFaceColor",'black',"MarkerSize",3)
hM_11=plot(real(M_eig(1,1)),imag(M_eig(1,1)),'o',"MarkerEdgeColor",[8,104,172]/255,"MarkerFaceColor",[8,104,172]/255,"MarkerSize",3);
hM_12=plot(real(M_eig(2,1)),imag(M_eig(2,1)),'o',"MarkerEdgeColor",[43,140,190]/255,"MarkerFaceColor",[43,140,190]/255,"MarkerSize",3);
hM_13=plot(real(M_eig(3,1)),imag(M_eig(3,1)),'o',"MarkerEdgeColor",'r',"MarkerFaceColor",'r',"MarkerSize",3);


xlabel('Real','fontsize',9);
ylabel('Imaginary','FontSize',9);
set(gca,'YLim',[-1.1 1.1]);
set(gca,'XLim',[-1.57 1.57]);
box on;

%%
subplot(1,2,2)
set(gca,'unit', 'centimeters', 'position', [8, 1, 5, 3.5]);
set(gca,'FontSize',8) %坐标轴字体的大小

%画不稳定的区域
point1=[-3 3 3 -3];
point2=[-3 -3 3 3];
Color=[0.75 0.75 0.75];
a=patch(point1,point2,Color) ;
alpha(0.6);
set(a,'edgecolor','none');
hold on;

% %画单位圆
% t=0:0.01:2*pi;
% x=cos(t);
% y=sin(t);
% h_circle=plot(x,y,'linestyle','-',"Color",'black','LineWidth',1.1)
% axis equal;

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
%colorbar 



%画不同半径的稳定性
% r=0.8;
r1=0.8;
theta=pi*(-1:0.01:1);
x_1=r1*cos(theta);
y_1=r1*sin(theta);
h_1=plot(x_1,y_1,'linestyle','--',"Color",[82 82 82]/255,'LineWidth',1)
stability_1=-log(r1)

%r=0.6;
r2=0.6;
theta=pi*(-1:0.01:1);
x_2=r2*cos(theta);
y_2=r2*sin(theta);
h_2=plot(x_2,y_2,'linestyle',':',"Color",[82 82 82]/255,'LineWidth',1)
stability_2=-log(r2)

%r=0.4;
r3=0.4;
theta=pi*(-1:0.01:1);
x_3=r3*cos(theta);
y_3=r3*sin(theta);
h_3=plot(x_3,y_3,'linestyle','-',"Color",[82 82 82]/255,'LineWidth',1)
stability_3=-log(r3)



%画群落矩阵的特征值分布
A=[-1,-0.6,0.3;
    0.6,-1,0.8;
    0.3,0.6,-1];

M=eye(3)+A;

M_eig=eig(M)

%hM_2=plot(real(M_eig),imag(M_eig),'o',"MarkerEdgeColor",'black',"MarkerFaceColor",'black',"MarkerSize",3)
hM_11=plot(real(M_eig(1,1)),imag(M_eig(1,1)),'o',"MarkerEdgeColor",[8,104,172]/255,"MarkerFaceColor",[8,104,172]/255,"MarkerSize",3);
hM_12=plot(real(M_eig(2,1)),imag(M_eig(2,1)),'o',"MarkerEdgeColor",[43,140,190]/255,"MarkerFaceColor",[43,140,190]/255,"MarkerSize",3);
hM_13=plot(real(M_eig(3,1)),imag(M_eig(3,1)),'o',"MarkerEdgeColor",'r',"MarkerFaceColor",'r',"MarkerSize",3);

legend([h_1 h_2 h_3],{ '0.22','0.51','0.92'},'fontsize',8)
legend boxoff


xlabel('Real','fontsize',9);
ylabel('Imaginary','FontSize',9);
set(gca,'YLim',[-1.1 1.1]);
set(gca,'XLim',[-1.57 1.57]);
box on;