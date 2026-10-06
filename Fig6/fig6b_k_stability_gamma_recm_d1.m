
% This file can be used to generate Fig. 6b 
% d=1




%% 横坐标为k 纵坐标为stability 随机、捕食-被捕食、竞争、互惠四条线画在一个小图里

%
clear all
clc



figure
set(gcf,'unit','centimeters','position',[15 5 10 6]) %长9 宽6
subplot(1,1,1)
set(gca,'unit', 'centimeters', 'position', [2 1 7 4]);
set(gca,'FontSize',8) %坐标轴字体的大小



%
point1=[0 3.2 3.2 0];
point2=[-2 -2 0 0];
Color=[0.75 0.75 0.75];
aa=patch(point1,point2,Color) 
hold on
alpha(0.6)
set(aa,'edgecolor','none')
text(0.1,-0.05,'Unstable','Color','black','FontSize',8);


d=1;
sigma=0.05;
S=100;  
C=0.2;
n=50;
%val_kesti=0.1:0.1:3.2; 
%val_ksimu=0.1:0.4:2.9;
val_kesti=0.1:0.1:3.1; 
val_ksimu=0.2:0.4:3;

%% random的三点估计法 
% i=1;
% for k=val_kesti
%     x_lm=-(d-1)-sigma*sqrt(S*C); %最左边特征值的实部
%     y_lm=0;   %最左边特征值的虚部
%     x_rm=-(d-1)+sigma*sqrt(S*C); %最右边特征值的实部
%     y_rm=0; %最右边特征值的虚部
%     x_um=-(d-1); %最上边特征值的实部
%     y_um=sigma*sqrt(S*C); %最上边特征值的虚部
% 
%     abs_lm=sqrt(x_lm^2+y_lm^2);%最左边特征值对应的绝对值
%     abs_rm=sqrt(x_rm^2+y_rm^2);%最右边特征值对应的绝对值
%     abs_um=sqrt(x_um^2+y_um^2);%最上边特征值对应的绝对值
%     J_abs=[abs_lm,abs_rm,abs_um]; %为了用max 构造矩阵 详情见max用法 
%     alpha_r(1,i)=-log(max(J_abs));
%     i=i+1;
% end

% Random 群落的精确理论值
i = 1;

for k = val_kesti

    % Gamma分布绝对值的均值
    E_abs = sigma*sqrt(k/(k+1));

    % --------------------------------------------------------
    % Random community
    % --------------------------------------------------------

    % 非对角元均值
    E = 0;

    % 非对角元方差
    V = C*sigma^2 - E^2;

    % M_ij 和 M_ji 的乘积均值
    R = 0;

    % 相关系数
    r = (R-E^2)/V;

    % --------------------------------------------------------
    % 椭圆半轴
    % --------------------------------------------------------

    a = sqrt(S*V)*(1+r);
    b = sqrt(S*V)*(1-r);

    % --------------------------------------------------------
    % 椭圆中心
    % --------------------------------------------------------

    c = 1-d-E;


    % --------------------------------------------------------
    % 椭圆的精确最大模
    % --------------------------------------------------------

    if a >= b

        % 第一种情况
        rho_ell = abs(c) + a;

    else

        % 此时 b > a
        condition = abs(a*c/(b^2-a^2));

        if condition <= 1

            % 第二种情况：最大模位于椭圆弧上
            rho_ell = b * sqrt( ...
                (b^2 + c^2 - a^2)/(b^2-a^2) ...
                );

        else

            % 第三种情况：最大模位于实轴端点
            rho_ell = abs(c) + a;

        end

    end


    % --------------------------------------------------------
    % 实离群特征值的模
    % --------------------------------------------------------

    rho_outlier = abs(1-d+(S-1)*E);


    % --------------------------------------------------------
    % 整个谱的最大模
    % --------------------------------------------------------

    rho_M = max(rho_ell,rho_outlier);


    % --------------------------------------------------------
    % alpha = -log(rho)
    % --------------------------------------------------------

    alpha_r(i) = -log(rho_M);

    i = i+1;

end

h_r=plot(val_kesti,alpha_r,'linestyle','-',"Color",'[0 0.4470 0.7410]','LineWidth',1)
hold on 

%random 数值仿真
m=1;
for k=val_ksimu
a=k;
b=sigma*sqrt(1/(k^2+k));
J=zeros(1,n);
for t=1:n
    M=zeros(S);
for i=1:S
    for j=1:S
        if i==j
            M(i,j)=-(d-1);
        else 
             p1=rand;
            if p1<=C
                p2=rand;
                if p2<0.5
                    M(i,j)=gamrnd(a,b);
                else
                    M(i,j)=-gamrnd(a,b);
                end
           end
        end
    end
end
J(1,t)=-log(max(abs(eig(M))));
end
J_mean(1,m)=mean(J); %random占J_mean的第一行
m=m+1;
end
%
h_rr=plot(val_ksimu,J_mean(1,:),'d','Markersize',4,'MarkerEdgeColor',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410])




hold on;

%% 捕食被捕食
% i=1;
% for k=val_kesti
% E_abs=sigma*sqrt(k/(k+1));
% x_lm=-(d-1)-sigma*sqrt(S*C)*(1-(E_abs^2)/(sigma^2));
% y_lm=0;
% x_rm=-(d-1)+sigma*sqrt(S*C)*(1-(E_abs^2)/(sigma^2));
% y_rm=0;
% x_um=-(d-1);
% y_um=sigma*sqrt(S*C)*(1+(E_abs^2)/(sigma^2));
% 
% abs_lm=sqrt(x_lm^2+y_lm^2);
% abs_rm=sqrt(x_rm^2+y_rm^2);
% abs_um=sqrt(x_um^2+y_um^2);
% 
% J_abs=[abs_lm,abs_rm,abs_um]; %为了用max 构造矩阵 详情见max用法 
% alpha_e(1,i)=-log(max(J_abs));
% i=i+1;
% end

% Exploitative 群落的精确理论值
i = 1;

for k = val_kesti

    % Gamma分布绝对值的均值
    E_abs = sigma*sqrt(k/(k+1));

    % --------------------------------------------------------
    % Exploitative community
    % --------------------------------------------------------

    % 非对角元均值
    E = 0;

    % 非对角元方差
    V = C*sigma^2 - E^2;

    % M_ij 和 M_ji 的乘积均值
    R = -C*E_abs^2;

    % 相关系数
    r = (R-E^2)/V;

    % --------------------------------------------------------
    % 椭圆半轴
    % --------------------------------------------------------

    a = sqrt(S*V)*(1+r);
    b = sqrt(S*V)*(1-r);

    % --------------------------------------------------------
    % 椭圆中心
    % --------------------------------------------------------

    c = 1-d-E;


    % --------------------------------------------------------
    % 椭圆的精确最大模
    % --------------------------------------------------------

    if a >= b

        % 第一种情况
        rho_ell = abs(c) + a;

    else

        % 此时 b > a
        condition = abs(a*c/(b^2-a^2));

        if condition <= 1

            % 第二种情况：最大模位于椭圆弧上
            rho_ell = b * sqrt( ...
                (b^2 + c^2 - a^2)/(b^2-a^2) ...
                );

        else

            % 第三种情况：最大模位于实轴端点
            rho_ell = abs(c) + a;

        end

    end


    % --------------------------------------------------------
    % 实离群特征值的模
    % --------------------------------------------------------

    rho_outlier = abs(1-d+(S-1)*E);


    % --------------------------------------------------------
    % 整个谱的最大模
    % --------------------------------------------------------

    rho_M = max(rho_ell,rho_outlier);


    % --------------------------------------------------------
    % alpha = -log(rho)
    % --------------------------------------------------------

    alpha_e(i) = -log(rho_M);

    i = i+1;

end

h_e=plot(val_kesti,alpha_e,'linestyle','-',"Color",[0.8500 0.3250 0.0980],'LineWidth',1)

%捕食被捕食的数值仿真
m=1;
for k=val_ksimu
a=k;
b=sigma*sqrt(1/(k^2+k));
J=zeros(1,n);
for t=1:n
    M=zeros(S);
for i=1:S
    for j=1:S
        if i==j
            M(i,j)=1-d;
        elseif i>j
                p1=rand;
            if p1<=C
                p2=rand;
                if p2<=0.5
                    M(i,j)=gamrnd(a,b);
                    M(j,i)=-gamrnd(a,b);
                else 
                    M(i,j)=-gamrnd(a,b);
                    M(j,i)=gamrnd(a,b);
                end
            end
        end
    end
end
J(1,t)=-log(max(abs(eig(M))));
end
J_mean(2,m)=mean(J);
m=m+1;
end
%
h_ee=plot(val_ksimu,J_mean(2,:),'o','Markersize',4,'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980])


%% 竞争的三点估计
% i=1;
% for k=val_kesti
%     
% E_abs=sigma*sqrt(k/(k+1));
% x_lm=-(d-1)-(S-1)*C*E_abs; 
% y_lm=0;
% x_rm=-(d-1)+C*E_abs+sqrt(S*C*(sigma^2-C*E_abs^2))*(sigma^2+(1-2*C)*E_abs^2)/(sigma^2-C*E_abs^2);
% y_rm=0;
% x_um=-(d-1)+C*E_abs;
% y_um=sqrt(S*C*(sigma^2-C*E_abs^2))*(sigma^2-E_abs^2)/(sigma^2-C*E_abs^2);
% 
% abs_lm=sqrt(x_lm^2+y_lm^2);
% abs_rm=sqrt(x_rm^2+y_rm^2);
% abs_um=sqrt(x_um^2+y_um^2);
% 
% J_abs=[abs_lm,abs_rm,abs_um];
% alpha_c(1,i)=-log(max(J_abs));
% i=i+1;
% end

% Competitive 群落的精确理论值
i = 1;

for k = val_kesti

    % Gamma分布绝对值的均值
    E_abs = sigma*sqrt(k/(k+1));

    % --------------------------------------------------------
    % Competitive community
    % --------------------------------------------------------

    % 非对角元均值
    E = -C*E_abs;

    % 非对角元方差
    V = C*sigma^2 - E^2;

    % M_ij 和 M_ji 的乘积均值
    R = C*E_abs^2;

    % 相关系数
    r = (R-E^2)/V;

    % --------------------------------------------------------
    % 椭圆半轴
    % --------------------------------------------------------

    a = sqrt(S*V)*(1+r);
    b = sqrt(S*V)*(1-r);

    % --------------------------------------------------------
    % 椭圆中心
    % --------------------------------------------------------

    c = 1-d-E;


    % --------------------------------------------------------
    % 椭圆的精确最大模
    % --------------------------------------------------------

    if a >= b

        % 第一种情况
        rho_ell = abs(c) + a;

    else

        % 此时 b > a
        condition = abs(a*c/(b^2-a^2));

        if condition <= 1

            % 第二种情况：最大模位于椭圆弧上
            rho_ell = b * sqrt( ...
                (b^2 + c^2 - a^2)/(b^2-a^2) ...
                );

        else

            % 第三种情况：最大模位于实轴端点
            rho_ell = abs(c) + a;

        end

    end


    % --------------------------------------------------------
    % 实离群特征值的模
    % --------------------------------------------------------

    rho_outlier = abs(1-d+(S-1)*E);


    % --------------------------------------------------------
    % 整个谱的最大模
    % --------------------------------------------------------

    rho_M = max(rho_ell,rho_outlier);


    % --------------------------------------------------------
    % alpha = -log(rho)
    % --------------------------------------------------------

    alpha_c(i) = -log(rho_M);

    i = i+1;

end

h_c=plot(val_kesti,alpha_c,'linestyle','-',"Color",[0.9290 0.6940 0.1250],'LineWidth',1)

%竞争的数值仿真
m=1;
for k=val_ksimu
a=k;
b=sigma*sqrt(1/(k^2+k));
J=zeros(1,n);
for t=1:n
    M=zeros(S);
for i=1:S
    for j=1:S
        if i==j
            M(i,j)=1-d;
        elseif i>j
                p=rand;
            if p<=C
                    M(i,j)=-gamrnd(a,b);
                    M(j,i)=-gamrnd(a,b);
            end
        end
    end
end
J(1,t)=-log(max(abs(eig(M))));
end
J_mean(3,m)=mean(J);
m=m+1;
end
%
h_cc=plot(val_ksimu,J_mean(3,:),'s','Markersize',5,'MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9290 0.6940 0.1250])




%% 互惠的三点估计
% i=1;
% for k=val_kesti
%     
% E_abs=sigma*sqrt(k/(k+1));
% x_lm=-(d-1)-C*E_abs-sqrt(S*C*(sigma^2-C*E_abs^2))*(sigma^2+(1-2*C)*E_abs^2)/(sigma^2-C*E_abs^2);
% y_lm=0;
% x_rm=-(d-1)+(S-1)*C*E_abs; 
% y_rm=0;
% x_um=-(d-1)-C*E_abs;
% y_um=sqrt(S*C*(sigma^2-C*E_abs^2))*(sigma^2-E_abs^2)/(sigma^2-C*E_abs^2);
% 
% abs_lm=sqrt(x_lm^2+y_lm^2);
% abs_rm=sqrt(x_rm^2+y_rm^2);
% abs_um=sqrt(x_um^2+y_um^2);
% 
% J_abs=[abs_lm,abs_rm,abs_um]; %为了用max 构造矩阵 详情见max用法 
% alpha_m(1,i)=-log(max(J_abs));
% i=i+1;
% end

% Mutualistic 群落的精确理论值
i = 1;

for k = val_kesti

    % Gamma分布绝对值的均值
    E_abs = sigma*sqrt(k/(k+1));

    % --------------------------------------------------------
    % Mutualistic community
    % --------------------------------------------------------

    % 非对角元均值
    E = C*E_abs;

    % 非对角元方差
    V = C*sigma^2 - E^2;

    % M_ij 和 M_ji 的乘积均值
    R = C*E_abs^2;

    % 相关系数
    r = (R-E^2)/V;

    % --------------------------------------------------------
    % 椭圆半轴
    % --------------------------------------------------------

    a = sqrt(S*V)*(1+r);
    b = sqrt(S*V)*(1-r);

    % --------------------------------------------------------
    % 椭圆中心
    % --------------------------------------------------------

    c = 1-d-E;


    % --------------------------------------------------------
    % 椭圆的精确最大模
    % --------------------------------------------------------

    if a >= b

        % 第一种情况
        rho_ell = abs(c) + a;

    else

        % 此时 b > a
        condition = abs(a*c/(b^2-a^2));

        if condition <= 1

            % 第二种情况：最大模位于椭圆弧上
            rho_ell = b * sqrt( ...
                (b^2 + c^2 - a^2)/(b^2-a^2) ...
                );

        else

            % 第三种情况：最大模位于实轴端点
            rho_ell = abs(c) + a;

        end

    end


    % --------------------------------------------------------
    % 实离群特征值的模
    % --------------------------------------------------------

    rho_outlier = abs(1-d+(S-1)*E);


    % --------------------------------------------------------
    % 整个谱的最大模
    % --------------------------------------------------------

    rho_M = max(rho_ell,rho_outlier);


    % --------------------------------------------------------
    % alpha = -log(rho)
    % --------------------------------------------------------

    alpha_m(i) = -log(rho_M);

    i = i+1;

end

h_m=plot(val_kesti,alpha_m,'linestyle','--',"Color",[0.4940 0.1840 0.5560],'LineWidth',1)

%互惠数值仿真
m=1;
for k=val_ksimu
a=k;
b=sigma*sqrt(1/(k^2+k));

J=zeros(1,n);
for t=1:n
    M=zeros(S);
for i=1:S
    for j=1:S
        if i==j
            M(i,j)=1-d;
        elseif i>j
                p=rand;
            if p<=C
                    M(i,j)=gamrnd(a,b);
                    M(j,i)=gamrnd(a,b);
            end
        end
    end
end
J(1,t)=-log(max(abs(eig(M))));
end
J_mean(4,m)=mean(J);
m=m+1;
end
%
k=val_ksimu;
h_mm=plot(val_ksimu,J_mean(4,:),'+','Color',[0.4940 0.1840 0.5560],'Markersize',4)
%
axis ([0,3.2,-0.15,1.65])

xlabel('$$k$$','Interpreter','latex','fontsize',9);
ylabel('Stability','FontSize',9); 
box on

% legend([h_r h_e h_c h_m h_rr h_ee h_cc h_mm],{'Random(Esti.)', 'Exploitative','Competitive','Mutualistic','Random', '+/-','-/-','+/+'},'FontSize',8)
% legend boxoff;%图例去边框

title('$$d=1$$', ...
    'Interpreter','latex', ...
    'FontSize',9);
