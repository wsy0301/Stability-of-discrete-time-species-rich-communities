

% This file can be used to generate Fig. 5a

clear all
clc

figure
set(gcf,'unit','centimeters','position',[15 5 12 9]) %长11 宽7
subplot(1,1,1)
set(gca,'unit', 'centimeters', 'position', [2 1 7 4]);
set(gca,'FontSize',8) %坐标轴字体的大小


%% 画不稳定的灰色区域
point1=[0 3 3 0];
point2=[-0.4 -0.4 0 0];
Color=[0.75 0.75 0.75];
aa=patch(point1,point2,Color) 

alpha(0.6)
set(aa,'edgecolor','none')

text(0.25,-0.21,'Unstable','Color','black','FontSize',8)
hold on 
%% 系统参数 随着自调节d的变化 稳定性指标的变化 稳定性指标为log10(max(svd(J_abs)))

sigma=0.05;
val_d=0:0.01:3;
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
E=0;
V=C*sigma^2;
R=0;%这里的R就是相关性 论文里的rho

a=sqrt(S*V)*(1+((R-E^2)/V));
b=sqrt(S*V)*(1-((R-E^2)/V));


i=1;
for d=val_d
    c=1-d-E;
lambda_r=abs(c)+a;
alpha_r(i)=-log(lambda_r); %求三个特征值对应的绝对值的最大——alpha
i=i+1;
end
val_d=0:0.01:3; 
h_r=plot(val_d,alpha_r,'linestyle','-',"Color",'[0 0.4470 0.7410]','LineWidth',1)
hold on 

%% 捕食-被捕食
E_abs=sigma*sqrt(2/pi);
E=0;
V=C*sigma^2-E^2;
R=-C*E_abs^2;%这里的R就是相关性 论文里的rho

a=sqrt(S*V)*(1+((R-E^2)/V));
b=sqrt(S*V)*(1-((R-E^2)/V));



i = 1;
for d=val_d
    
    % 椭圆中心
    c=1-d-E;

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
            % 第二种情况
            rho_ell = b * sqrt((b^2 + c^2 - a^2)/(b^2-a^2));
        else
            % 第三种情况
            rho_ell = abs(c) + a;
        end
    end

    % --------------------------------------------------------
    % 实离群特征值的模
    % --------------------------------------------------------
    rho_outlier = abs(1 - d + (S-1)*E);
    % --------------------------------------------------------
    % 整个谱的最大模
    % --------------------------------------------------------
    rho_M = max(rho_ell, rho_outlier);    
    % --------------------------------------------------------
    % alpha = -log(rho)
    % --------------------------------------------------------
    alpha_e(i) = -log(rho_M);
    i = i + 1;
end

val_d=0:0.01:3; 
h_e=plot(val_d,alpha_e,'linestyle','-',"Color",'[0.8500 0.3250 0.0980]','LineWidth',1)


%% 独立竞争

E_abs=sigma*sqrt(2/pi);
E=-C*E_abs;
V=C*sigma^2-E^2;
R=C*E_abs^2;%这里的R就是相关性 论文里的rho

a=sqrt(S*V)*(1+((R-E^2)/V));
b=sqrt(S*V)*(1-((R-E^2)/V));



i = 1;
for d=val_d
    
    % 椭圆中心
    c=1-d-E;


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
            % 第二种情况
            rho_ell = b * sqrt((b^2 + c^2 - a^2)/(b^2-a^2));
        else
            % 第三种情况
            rho_ell = abs(c) + a;
        end
    end

    % --------------------------------------------------------
    % 实离群特征值的模
    % --------------------------------------------------------
    rho_outlier = abs(1 - d + (S-1)*E);
    % --------------------------------------------------------
    % 整个谱的最大模
    % --------------------------------------------------------
    rho_M = max(rho_ell, rho_outlier);
    % --------------------------------------------------------
    % alpha = -log(rho)
    % --------------------------------------------------------
    alpha_c(i) = -log(rho_M);
    i = i + 1;
end

val_d=0:0.01:3; 
h_c=plot(val_d,alpha_c,'linestyle','-',"Color",'[0.9290 0.6940 0.1250]','LineWidth',1)


%% 独立互惠 +/+
E_abs=sigma*sqrt(2/pi);
E=C*E_abs;
V=C*sigma^2-E^2;
R=C*E_abs^2;%这里的R就是相关性 论文里的rho

a=sqrt(S*V)*(1+((R-E^2)/V));
b=sqrt(S*V)*(1-((R-E^2)/V));
c=1-d-E;


i = 1;
for d=val_d
    
    % 椭圆中心
    c=1-d-E;


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
            % 第二种情况
            rho_ell = b * sqrt((b^2 + c^2 - a^2)/(b^2-a^2));
        else
            % 第三种情况
            rho_ell = abs(c) + a;
        end
    end

    % --------------------------------------------------------
    % 实离群特征值的模
    % --------------------------------------------------------
    rho_outlier = abs(1 - d + (S-1)*E);
    % --------------------------------------------------------
    % 整个谱的最大模
    % --------------------------------------------------------
    rho_M = max(rho_ell, rho_outlier);
    % --------------------------------------------------------
    % alpha = -log(rho)
    % --------------------------------------------------------
    alpha_m(i) = -log(rho_M);
    i = i + 1;
end


val_d=0:0.01:3; 
h_m=plot(val_d,alpha_m,'linestyle','-',"Color",'[0.4940 0.1840 0.5560]','LineWidth',1)


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


%% 生成50个捕食-被捕食者矩阵
% 每一个矩阵求特征值——求特征值的模--取最大——50个求均值

% n=50;
% J=cell(n,1);
% M=zeros(S,S);
% for d=0.1:0.1:3;
% for t=1:n
% for i=1:S
%     for j=1:S
%         if i==j
%             M(i,j)=-(d-1);
%         else
%             if i>j
%                 p1=rand;
%             if p1<=C
%                 p2=rand;
%                 if p2<=0.5
%                     M(i,j)=abs(normrnd(0,sigma,1,1));
%                     M(j,i)=-abs(normrnd(0,sigma,1,1));
%                 else 
%                     M(i,j)=-abs(normrnd(0,sigma,1,1));
%                     M(j,i)=abs(normrnd(0,sigma,1,1));
%                 end
%             end
%             end
%         end
%     end
% end
% J{t,1}=M;
% end
% val_z=cell(S,1);
% alpha1=cell(n,1);
% for t=1:n
%     lam=eig(cell2mat(J{t,1}));%求群落矩阵M的特征值 
% for k=1:S
%     val_z{1,k}=abs(lam(k)); %求M矩阵的特征值的模
% end
% alpha1{t,1}=-log(max(cell2mat(val_z)));%每一个群落矩阵M中S个特征值对应的特征根的最大实部——取复数——alpha
% end
% alpha_mean=mean(cell2mat(alpha1));
% 
% h_ee=plot(d,alpha_mean,'o','Markersize',4,'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980])
% end


n = 50;                       % 每个 d 下随机群落数
d_all = 0.1:0.1:3;            % d 的取值

alpha_mean = zeros(size(d_all));

for q = 1:length(d_all)

    d = d_all(q);

    alpha1 = zeros(n,1);

    for t = 1:n

        % =====================================================
        % 每生成一个新的群落矩阵，都必须重新初始化 M
        % =====================================================
        M = zeros(S,S);

        % 对角线
        for i = 1:S
            M(i,i) = -(d-1);
        end

        % =====================================================
        % 生成 exploitative interaction
        % =====================================================
        for i = 2:S
            for j = 1:i-1

                p1 = rand;

                if p1 <= C

                    p2 = rand;

                    if p2 <= 0.5

                        M(i,j) = abs(normrnd(0,sigma));
                        M(j,i) = -abs(normrnd(0,sigma));

                    else

                        M(i,j) = -abs(normrnd(0,sigma));
                        M(j,i) = abs(normrnd(0,sigma));

                    end

                end

            end
        end


        % =====================================================
        % 特征值
        % =====================================================
        lam = eig(M);

        % 最大模特征值
        rho_M = max(abs(lam));

        % decay rate alpha
        alpha1(t) = -log(rho_M);

    end


    % =========================================================
    % 当前 d 下 50 个群落的平均 alpha
    % =========================================================
    alpha_mean(q) = mean(alpha1);

end


% =============================================================
% 最后统一作图
% =============================================================
h_ee = plot(d_all,alpha_mean,'o', ...
    'MarkerSize',4, ...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], ...
    'MarkerFaceColor',[0.8500 0.3250 0.0980]);

hold on

%% 竞争 50个矩阵 每一个矩阵求特征值——求特征值的模--取最大——50个求均值 


% n=50;
% J=cell(n,1);
% 
% for d=0.1:0.1:3;
% for t=1:n
% for i=1:S
%     for j=1:S
%         if i==j
%             M{i,j}=-(d-1);
%         else if i>j
%                 p=rand;
%             if p<=C
%                     M{i,j}=-abs(normrnd(0,sigma,1,1));
%                     M{j,i}=-abs(normrnd(0,sigma,1,1));
%                 else 
%                     M{j,i}=0;
%                     M{i,j}=0;
%             end
%             end
%         end
%     end
% end
% J{t,1}=M;
% end
% val_z=cell(S,1);
% alpha1=cell(n,1);
% for t=1:n
%     lam=eig(cell2mat(J{t,1}));%求群落矩阵M的特征值 
% for k=1:S
%     val_z{1,k}=abs(lam(k)); %求M矩阵的特征值的模
% end
% alpha1{t,1}=-log(max(cell2mat(val_z)));%每一个群落矩阵M中S个特征值对应的特征根的最大实部——取复数——alpha
% end
% alpha_mean=mean(cell2mat(alpha1));
% 
% h_cc=plot(d,alpha_mean,'s','Markersize',5,'MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9290 0.6940 0.1250])
% end

n = 50;                       % 每个 d 下随机群落数
d_all = 0.1:0.1:3;            % d 的取值

alpha_mean_c = zeros(size(d_all));

for q = 1:length(d_all)

    d = d_all(q);

    alpha1 = zeros(n,1);

    for t = 1:n

        % =====================================================
        % 每生成一个新的群落矩阵，都重新初始化 M
        % =====================================================
        M = zeros(S,S);

        % 对角线
        for i = 1:S
            M(i,i) = -(d-1);
        end


        % =====================================================
        % 生成 competitive interaction
        % 两个方向均为负
        % =====================================================
        for i = 2:S
            for j = 1:i-1

                p = rand;

                if p <= C

                    M(i,j) = -abs(normrnd(0,sigma));
                    M(j,i) = -abs(normrnd(0,sigma));

                end

            end
        end


        % =====================================================
        % 特征值
        % =====================================================
        lam = eig(M);

        % 最大模特征值，即 spectral radius
        rho_M = max(abs(lam));

        % decay rate alpha
        alpha1(t) = -log(rho_M);

    end


    % =========================================================
    % 当前 d 下 50 个群落的平均 alpha
    % =========================================================
    alpha_mean_c(q) = mean(alpha1);

end


% =============================================================
% 最后统一作图
% =============================================================
h_cc = plot(d_all,alpha_mean_c,'s', ...
    'MarkerSize',5, ...
    'MarkerEdgeColor',[0.9290 0.6940 0.1250], ...
    'MarkerFaceColor',[0.9290 0.6940 0.1250]);

hold on

%% 互惠 50个矩阵 每一个矩阵求特征值——求对应的特征根——最大实部——50个求均值 

% n=50;
% J=cell(n,1);
% 
% for d=0.1:0.1:3;
% for t=1:n
% for i=1:S
%     for j=1:S
%         if i==j
%             M{i,j}=-(d-1);
%         else if i>j
%                 p=rand;
%             if p<=C
%                     M{i,j}=abs(normrnd(0,sigma,1,1));
%                     M{j,i}=abs(normrnd(0,sigma,1,1));
%                 else 
%                     M{j,i}=0;
%                     M{i,j}=0;
%             end
%             end
%         end
%     end
% end
% J{t,1}=M;
% end
% val_z=cell(S,1);
% alpha1=cell(n,1);
% for t=1:n
%     lam=eig(cell2mat(J{t,1}));%求群落矩阵M的特征值 
% for k=1:S
%    val_z{1,k}=abs(lam(k)); %求M矩阵的特征值的模
% end
% alpha1{t,1}=-log(max(cell2mat(val_z)));%每一个群落矩阵M中S个特征值对应的特征根的最大实部——取复数——alpha
% end
% alpha_mean=mean(cell2mat(alpha1));
% 
% h_mm=plot(d,alpha_mean,'+','Color',[0.4940 0.1840 0.5560],'Markersize',4)
% end

n = 50;                       % 每个 d 下随机群落数
d_all = 0.1:0.1:3;            % d 的取值

alpha_mean_m = zeros(size(d_all));

for q = 1:length(d_all)

    d = d_all(q);

    alpha1 = zeros(n,1);

    for t = 1:n

        % =====================================================
        % 每生成一个新的群落矩阵，都重新初始化 M
        % =====================================================
        M = zeros(S,S);

        % 对角线
        for i = 1:S
            M(i,i) = -(d-1);
        end


        % =====================================================
        % 生成 mutualistic interaction
        % 两个方向均为正
        % =====================================================
        for i = 2:S
            for j = 1:i-1

                p = rand;

                if p <= C

                    M(i,j) = abs(normrnd(0,sigma));
                    M(j,i) = abs(normrnd(0,sigma));

                end

            end
        end


        % =====================================================
        % 特征值
        % =====================================================
        lam = eig(M);

        % 最大模特征值，即 spectral radius
        rho_M = max(abs(lam));

        % decay rate alpha
        alpha1(t) = -log(rho_M);

    end


    % =========================================================
    % 当前 d 下 50 个群落的平均 alpha
    % =========================================================
    alpha_mean_m(q) = mean(alpha1);

end


% =============================================================
% 最后统一作图
% =============================================================
h_mm = plot(d_all,alpha_mean_m,'+', ...
    'Color',[0.4940 0.1840 0.5560], ...
    'MarkerSize',4);

hold on


%% 
axis ([0,2.2,-0.3,2.3])

xlabel('$$s$$','Interpreter','latex','fontsize',9);
ylabel('Stability','FontSize',9); 
box on


%%
%legend([h_r h_e h_c h_m h_rr h_ee h_cc h_mm],{'Random(Esti.)', '+/-','-/-','+/+','Random', '+/-','-/-','+/+'},'FontSize',8)
%legend boxoff;%图例去边框
%%
%title('sigma=0.05,S=100,C=0.1') % 标题
%title('sigma=0.05,S=100,C=0.2') % 标题
%title('sigma=0.05,S=100,C=0.3') % 标题
%title('sigma=0.05,S=200,C=0.1') % 标题
%title('sigma=0.05,S=300,C=0.1') % 标题

%legend([h_r h_e h_mix h_c h_m h_rr h_ee h_mixmix h_cc h_mm],{'Random(Esti.)', '+/-(Esti)','Mix(Esti.)','-/-(Esti.)','+/+(Esti.)','Random(Simu.)', '+/-(Simu.)','Mix(Simu.)','-/-(Simu.)','+/+(Simu.)'})