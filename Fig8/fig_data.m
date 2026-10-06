% This file can be used to generate Fig. 8

clear all
clc

figure 
set(gcf,'unit','centimeters','position',[5 5 30 9]) %举例左侧10  举例下侧5  

%%
%相互作用的参数
sigma=0.05;


%% 蝗虫-植物 蝗虫-蝗虫 正态分布 也是正文中的图 
subplot(1,4,1)
set(gca,'unit', 'centimeters', 'position', [2 1 5 4]);
set(gca,'FontSize',8) %坐标轴字体的大小


%%% 不稳定的区域 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 point1=[-4 4 4 -4];
 point2=[-2 -2 0 0];
 Color=[0.75 0.75 0.75];
 a=patch(point1,point2,Color) 
 alpha(0.6)
 set(a,'edgecolor','none')
 
 text(0.1,-0.3,'Unstable','Color','black','FontSize',8)
 hold on
 
 
% friedman,j gore j的17年nee

%%% 蝗虫和小麦 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B_connect1=[1	0	1;
% 1	0	1;
% 1	0	0;
% 0	1	0];

B_connect1=[1	0;
1	0;
1	0;
0	1];


%
[B_row,B_column]=size(B_connect1); %获取矩阵的行数和列数
B_connect1(B_connect1>0)=1;
B_connect2=B_connect1';
%
m=1;
for d=0:0.1:2.3
    
    M_triu=cell(B_row,B_column);
% M的上三角矩阵  每个元素都来自半正态分布
for i=1:B_row
    for j=1:B_column
     M_triu{i,j}=abs(normrnd(0,sigma,1,1)); %M的上三角矩阵 
     end
end
M_triu=cell2mat(M_triu).*B_connect1;

 M_tril=cell(B_column,B_row);
% M的下三角矩阵  每个元素都来自半正态分布
for i=1:B_column
    for j=1:B_row
     M_tril{i,j}=-abs(normrnd(0,sigma,1,1)); %M的下三角矩阵 
     end
end
M_tril=cell2mat(M_tril).*B_connect2;

%M的(1,1)块矩阵
M_11=zeros(B_row);
M_11(logical(eye(size(M_11))))=(1-d)*ones(B_row,1);

%M的(2,2)块矩阵
M_22=zeros(B_column);
M_22(logical(eye(size(M_22))))=(1-d)*ones(B_column,1);

M=[M_11,M_triu;M_tril,M_22];

J(1,m)=-log(max(abs(eig(M))));
m=m+1;
end
d=0:0.1:2.3;
h_e=plot(d,J,'o','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor','none','Markersize',4)
%h_e=scatter(d,J,5,'o','MarkerEdgeColor',[0.8500 0.3250 0.0980])
%h_e=plot(d,J,'-*',"Color",'r','LineWidth',1)
hold on



 %分段拟合
 %先拟合d=0到1 共11个点 
 p_e1=polyfit(d(1:11),J(1:11),3);
 plot(d(1:11),polyval(p_e1,d(1:11)),"Color",[0.8500 0.3250 0.0980],'LineWidth',1);
 %再拟合d=1之后的 共14个点 第112到第24
 p_e1=polyfit(d(11:24),J(11:24),3);
 plot(d(11:24),polyval(p_e1,d(11:24)),"Color",[0.8500 0.3250 0.0980],'LineWidth',1);
% 
%p_e=polyfit(d,J,6);
%plot(d,polyval(p_e,d),"Color",[0.8500 0.3250 0.0980],'LineWidth',1);

%%% 蝗虫和蝗虫之间的竞争 选取自altuda
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 先从植物-蝗虫的邻接网络中提取蝗虫-蝗虫之间的交互
P=[0	0	1	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0;
0	0	1	1	1	0	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	1	0	1	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	1	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
1	0	1	0	1	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	1	0	1	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	1	0	1	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0;
0	0	1	1	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0;
0	0	1	0	1	1	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	1	0	1	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	1	0	1	0	0	0	0	0	0	0	0	1	0	0	0	0	1	0	1	1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	1	1	0	0	1	1;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	1	0	0	0	0	0	0	0	0	0	1	0	0	0	1	0	0	0	1	1	1	0	0	1	0	0	0	0	0	0	0	0	0;
0	0	1	0	0	0	0	0	0	0	0	1	0	1	0	1	1	1	1	0	0	0	1	1	1	0	0	1	1	0	1	0	1	1	0	0	0	0	1	1	1	1	0	1	0	0	0	1	1	0	0	0;
0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	1	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	0	0	0	0	0	0	0	0	0	1	0	0	0	1	0	0	1	0	1	1	0	0	0	0	0	0	1	1	0	0	0	1	0	0	1	1	0	1	0	1	0	1	0	0	0	0	0	0;
0	0	1	0	0	0	0	0	0	0	0	1	0	1	1	0	0	1	1	0	0	0	1	1	0	1	0	0	0	0	0	0	0	0	0	1	0	1	1	0	1	0	1	0	0	0	0	0	0	0	0	0;
0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	1	1	0	0	0	0	1	1	0	0	1	0	0	0	0	1	0	1	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0;
0	1	1	0	0	1	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	1	1	1	0	0	0	0	0	0	0	0	0	0	0	1	0	1	1	0	0	0	1	0	1	0	0	0	0	0	0	0;
0	0	1	0	1	1	0	1	1	0	0	0	0	0	0	1	0	0	0	0	0	0	0	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0];


K = P * P';   

% 2. 清除自连接（对角线）
K(1:size(K,1)+1:end) = 0; 

% 3. 将有连接的地方设为1
K = double(K > 0);

[K_row,K_column]=size(K); %获取矩阵的行数和列数



m=1;
for d=0:0.1:2.3
    A=zeros(K_row,K_column);
    for i=1:K_row
    for j=1:K_column
     A(i,j)=-abs(normrnd(0,sigma,1,1)); 
    end
    end
    M=A.*K;
    M(logical(eye(size(M))))=(1-d)*ones(K_row,1);
J(1,m)=-log(max(abs(eig(M))));
m=m+1;
end

d=0:0.1:2.3;
h_c=plot(d,J,'s','MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerFaceColor','none','Markersize',4.5)
hold on


%分段拟合
%先拟合d=0到1 共11个点 
%p_m1=polyfit(d(1:11),J(1:11),3);
%plot(d(1:11),polyval(p_m1,d(1:11)),"Color",[0.9290 0.6940 0.1250],'LineWidth',1);
%再拟合d=1之后的 共14个点 第112到第24
%p_m1=polyfit(d(11:24),J(11:24),3);
%plot(d(11:24),polyval(p_m1,d(11:24)),"Color",[0.9290 0.6940 0.1250],'LineWidth',1);


 p_m=polyfit(d,J,6);
 plot(d,polyval(p_m,d),"Color",[0.9290 0.6940 0.1250],'LineWidth',1);



 axis ([-0.05,2.35,-0.8,3.5])
%  
 xlabel('$$s$$','Interpreter','latex','fontsize',9);
 ylabel('Stability','FontSize',9); 
 box on 
 legend([h_e,h_c],{'Grasshoppers-Plants (+/-)','Grasshoppers (-/-)'},'fontsize',8)
 legend boxoff;%图例去边框
 title('Normal distribution','FontSize',9);
 
 
 %% 先生成d较小时 真实群落的稳定性排序 d=0.5 
 %正文图4的真实数据
%我们为了说明random是最稳定的 


P=[0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0;
0	1	1	1	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0;
0	0	1	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	0	0	1	1	0	1	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	1	1	0	0	0	0	0	1	1	0	1	0;
0	0	1	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
1	0	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0;
0	0	1	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	1	0	0	1	0	1	0	0;
0	0	1	1	1	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	1	0	1	1	0	1	1	0	0	0	0	1	0	1	0	1	1	0	1	0	0	1	1	0	0	0	0	0	1	0	0	1	0;
0	1	1	1	1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0;
0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	0	1	1	0	0	0	1	1	1	0	0	0	0	0	1	1	0	1	0	1	0	0	1	1	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	1	1	0	0	0	1	1	1	0	0	0;
0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	1	1	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0;
0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	1	0	0	0	0	0	0	1	0	1	1	1	1	0	0	0	1	1	0	0	0;
0	0	0	1	0	0	0	0	0	0	0	0	1	1	0	0	1	0	0	1	0	1	0	0	0	1	0	1	0	1	1	1	0	1	0	1	0	1	1	1	1	0	1	1	0	0	0	0	0	1	1	1	0	1;
0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	1	0	0	0	1	0	0	0	0	0	1	0	0	0	1	0	0	0	0	1	1	0	0	0	0	0	0	1	0	0	0;
0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0];



A = P * P';   

% 2. 清除自连接（对角线）
A(1:size(A,1)+1:end) = 0; 

% 3. 将有连接的地方设为1
A = double(A > 0)

%根据M生成的连接矩阵K
K=A;

d=0.5;

%均值和方差
sigma=0.01;

%规模
S=size(K,1);

% 借助网络结构生成random的
for k=1:100
M_r=zeros(S);
for i=1:S
    for j=1:S
        M_r(i,j)=normrnd(0,sigma,1,1);
    end
end
M_r=M_r.*K;
M_r(logical(eye(size(M_r))))=(1-d)*ones(size(M_r,1),1);
stability_r(1,k)=-log(max(abs(eig(M_r))));
k=k+1;
end


      

% 捕食被捕食
for k=1:100
M_e=zeros(S);
for i=1:S
    for j=1:S
         if i>j
              p2=rand;
                if p2<=0.5
                    M_e(i,j)=abs(normrnd(0,sigma,1,1));
                    M_e(j,i)=-abs(normrnd(0,sigma,1,1));
                else 
                    M_e(i,j)=-abs(normrnd(0,sigma,1,1));
                    M_e(j,i)=abs(normrnd(0,sigma,1,1));
                end
         end
    end
end
M_e=M_e.*K;
M_e(logical(eye(size(M_e))))=(1-d)*ones(size(M_e,1),1);
stability_e(1,k)=-log(max(abs(eig(M_e))));
k=k+1;
end


% 竞争
for k=1:100
M_c=zeros(S);
for i=1:S
    for j=1:S
         if i>j
                 M_c(i,j)=-abs(normrnd(0,sigma,1,1));
                 M_c(j,i)=-abs(normrnd(0,sigma,1,1));
         end
    end
end
M_c=M_c.*K;
M_c(logical(eye(size(M_c))))=(1-d)*ones(size(M_c,1),1);
stability_c(1,k)=-log(max(abs(eig(M_c))));
k=k+1;
end

% 互惠
for k=1:100
M_m=zeros(S);
for i=1:S;
    for j=1:S;
         if i>j
                 M_m(i,j)=abs(normrnd(0,sigma,1,1));
                 M_m(j,i)=abs(normrnd(0,sigma,1,1));
         end
    end
end
M_m=M_m.*K;
M_m(logical(eye(size(M_m))))=(1-d)*ones(size(M_m,1),1);
stability_m(1,k)=-log(max(abs(eig(M_m))));
k=k+1;
end


% 画图


subplot(1,4,2)
set(gca,'unit', 'centimeters', 'position', [8, 1, 5, 4]);

set(gca,'FontSize',8)
hold on 
stability=[mean(stability_r) mean(stability_e) mean(stability_c) mean(stability_m)];

b=bar(stability,0.4)
b.FaceColor = 'flat'; % 将 FaceColor 属性设置为 'flat' 以启用逐柱子颜色设置
b.CData(1, :) = [0 0.4470 0.7410]; % 第一个柱子的颜色
b.CData(2, :) = [0.8500 0.3250 0.0980]; % 第二个柱子的颜色
b.CData(3, :) = [0.9290 0.6940 0.1250]; % 第三个柱子的颜色
b.CData(4, :) = [0.4940 0.1840 0.5560]; % 第四个柱子的颜色

% set(b(1),'FaceColor',[0 0.4470 0.7410]);
% set(b(2),'FaceColor',[0.8500 0.3250 0.0980]);
% set(b(3),'FaceColor',[0.9290 0.6940 0.1250]);
% set(b(4),'FaceColor',[0.4940 0.1840 0.5560]);
set(b,'edgecolor','none');
%alpha(0.85);
axis([0.5 4.5 0 0.95]);
set(gca, 'xticklabel',{'Random','+/-','-/-','+/+'},'FontSize',8);
box on

ylabel('Stability','FontSize',9); 
xlabel('Interaction type','fontsize',9);
%title('Empirical dataset','fontsize',9);
title('$$d=0.5$$','Interpreter','latex','fontsize',9)

%%  生成d适中时 真实群落的稳定性排序 d=1 
 %正文图4的真实数据
%我们为了说明random是最稳定的 


P=[0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0;
0	1	1	1	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0;
0	0	1	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	0	0	1	1	0	1	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	1	1	0	0	0	0	0	1	1	0	1	0;
0	0	1	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
1	0	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0;
0	0	1	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	1	0	0	1	0	1	0	0;
0	0	1	1	1	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	1	0	1	1	0	1	1	0	0	0	0	1	0	1	0	1	1	0	1	0	0	1	1	0	0	0	0	0	1	0	0	1	0;
0	1	1	1	1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0;
0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	0	1	1	0	0	0	1	1	1	0	0	0	0	0	1	1	0	1	0	1	0	0	1	1	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	1	1	0	0	0	1	1	1	0	0	0;
0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	1	1	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0;
0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	1	0	0	0	0	0	0	1	0	1	1	1	1	0	0	0	1	1	0	0	0;
0	0	0	1	0	0	0	0	0	0	0	0	1	1	0	0	1	0	0	1	0	1	0	0	0	1	0	1	0	1	1	1	0	1	0	1	0	1	1	1	1	0	1	1	0	0	0	0	0	1	1	1	0	1;
0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	1	0	0	0	1	0	0	0	0	0	1	0	0	0	1	0	0	0	0	1	1	0	0	0	0	0	0	1	0	0	0;
0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0];



A = P * P';   

% 2. 清除自连接（对角线）
A(1:size(A,1)+1:end) = 0; 

% 3. 将有连接的地方设为1
A = double(A > 0)

%根据M生成的连接矩阵K
K=A;

d=1;

%规模
S=size(K,1);

% 借助网络结构生成random的
for k=1:100
M_r=zeros(S);
for i=1:S
    for j=1:S
        M_r(i,j)=normrnd(0,sigma,1,1);
    end
end
M_r=M_r.*K;
M_r(logical(eye(size(M_r))))=(1-d)*ones(size(M_r,1),1);
stability_r(1,k)=-log(max(abs(eig(M_r))));
k=k+1;
end


      

% 捕食被捕食
for k=1:100
M_e=zeros(S);
for i=1:S
    for j=1:S
         if i>j
              p2=rand;
                if p2<=0.5
                    M_e(i,j)=abs(normrnd(0,sigma,1,1));
                    M_e(j,i)=-abs(normrnd(0,sigma,1,1));
                else 
                    M_e(i,j)=-abs(normrnd(0,sigma,1,1));
                    M_e(j,i)=abs(normrnd(0,sigma,1,1));
                end
         end
    end
end
M_e=M_e.*K;
M_e(logical(eye(size(M_e))))=(1-d)*ones(size(M_e,1),1);
stability_e(1,k)=-log(max(abs(eig(M_e))));
k=k+1;
end


% 竞争
for k=1:100
M_c=zeros(S);
for i=1:S
    for j=1:S
         if i>j
                 M_c(i,j)=-abs(normrnd(0,sigma,1,1));
                 M_c(j,i)=-abs(normrnd(0,sigma,1,1));
         end
    end
end
M_c=M_c.*K;
M_c(logical(eye(size(M_c))))=(1-d)*ones(size(M_c,1),1);
stability_c(1,k)=-log(max(abs(eig(M_c))));
k=k+1;
end

% 互惠
for k=1:100
M_m=zeros(S);
for i=1:S;
    for j=1:S;
         if i>j
                 M_m(i,j)=abs(normrnd(0,sigma,1,1));
                 M_m(j,i)=abs(normrnd(0,sigma,1,1));
         end
    end
end
M_m=M_m.*K;
M_m(logical(eye(size(M_m))))=(1-d)*ones(size(M_m,1),1);
stability_m(1,k)=-log(max(abs(eig(M_m))));
k=k+1;
end


% 画图


subplot(1,4,3)
set(gca,'unit', 'centimeters', 'position', [14, 1, 5, 4]);

set(gca,'FontSize',8)
hold on 
stability=[mean(stability_r) mean(stability_e) mean(stability_c) mean(stability_m)];

b=bar(stability,0.4)
b.FaceColor = 'flat'; % 将 FaceColor 属性设置为 'flat' 以启用逐柱子颜色设置
b.CData(1, :) = [0 0.4470 0.7410]; % 第一个柱子的颜色
b.CData(2, :) = [0.8500 0.3250 0.0980]; % 第二个柱子的颜色
b.CData(3, :) = [0.9290 0.6940 0.1250]; % 第三个柱子的颜色
b.CData(4, :) = [0.4940 0.1840 0.5560]; % 第四个柱子的颜色

% set(b(1),'FaceColor',[0 0.4470 0.7410]);
% set(b(2),'FaceColor',[0.8500 0.3250 0.0980]);
% set(b(3),'FaceColor',[0.9290 0.6940 0.1250]);
% set(b(4),'FaceColor',[0.4940 0.1840 0.5560]);
set(b,'edgecolor','none');
%alpha(0.85);
axis([0.5 4.5 0 3.75]);
set(gca, 'xticklabel',{'Random','+/-','-/-','+/+'},'FontSize',8);
box on

ylabel('Stability','FontSize',9); 
xlabel('Interaction type','fontsize',9);
%title('Empirical dataset','fontsize',9);
title('$$d=1$$','Interpreter','latex','fontsize',9)

%%  d较大时 真实群落的稳定性排序 d=1.5 
 %正文图4的真实数据
%我们为了说明random是最稳定的 


P=[0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0;
0	1	1	1	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0;
0	0	1	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	0	0	1	1	0	1	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	1	1	0	0	0	0	0	1	1	0	1	0;
0	0	1	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
1	0	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0;
0	0	1	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	1	0	0	1	0	1	0	0;
0	0	1	1	1	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	1	0	1	1	0	1	1	0	0	0	0	1	0	1	0	1	1	0	1	0	0	1	1	0	0	0	0	0	1	0	0	1	0;
0	1	1	1	1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0;
0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	0	1	1	0	0	0	1	1	1	0	0	0	0	0	1	1	0	1	0	1	0	0	1	1	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	1	1	0	0	0	1	1	1	0	0	0;
0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	1	1	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0;
0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	1	0	0	0	0	0	0	1	0	1	1	1	1	0	0	0	1	1	0	0	0;
0	0	0	1	0	0	0	0	0	0	0	0	1	1	0	0	1	0	0	1	0	1	0	0	0	1	0	1	0	1	1	1	0	1	0	1	0	1	1	1	1	0	1	1	0	0	0	0	0	1	1	1	0	1;
0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	1	0	0	0	1	0	0	0	0	0	1	0	0	0	1	0	0	0	0	1	1	0	0	0	0	0	0	1	0	0	0;
0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0];



A = P * P';   

% 2. 清除自连接（对角线）
A(1:size(A,1)+1:end) = 0; 

% 3. 将有连接的地方设为1
A = double(A > 0)

%根据M生成的连接矩阵K
K=A;

d=1.5;

%规模
S=size(K,1);

% 借助网络结构生成random的
for k=1:100
M_r=zeros(S);
for i=1:S
    for j=1:S
        M_r(i,j)=normrnd(0,sigma,1,1);
    end
end
M_r=M_r.*K;
M_r(logical(eye(size(M_r))))=(1-d)*ones(size(M_r,1),1);
stability_r(1,k)=-log(max(abs(eig(M_r))));
k=k+1;
end


      

% 捕食被捕食
for k=1:100
M_e=zeros(S);
for i=1:S
    for j=1:S
         if i>j
              p2=rand;
                if p2<=0.5
                    M_e(i,j)=abs(normrnd(0,sigma,1,1));
                    M_e(j,i)=-abs(normrnd(0,sigma,1,1));
                else 
                    M_e(i,j)=-abs(normrnd(0,sigma,1,1));
                    M_e(j,i)=abs(normrnd(0,sigma,1,1));
                end
         end
    end
end
M_e=M_e.*K;
M_e(logical(eye(size(M_e))))=(1-d)*ones(size(M_e,1),1);
stability_e(1,k)=-log(max(abs(eig(M_e))));
k=k+1;
end


% 竞争
for k=1:100
M_c=zeros(S);
for i=1:S
    for j=1:S
         if i>j
                 M_c(i,j)=-abs(normrnd(0,sigma,1,1));
                 M_c(j,i)=-abs(normrnd(0,sigma,1,1));
         end
    end
end
M_c=M_c.*K;
M_c(logical(eye(size(M_c))))=(1-d)*ones(size(M_c,1),1);
stability_c(1,k)=-log(max(abs(eig(M_c))));
k=k+1;
end

% 互惠
for k=1:100
M_m=zeros(S);
for i=1:S;
    for j=1:S;
         if i>j
                 M_m(i,j)=abs(normrnd(0,sigma,1,1));
                 M_m(j,i)=abs(normrnd(0,sigma,1,1));
         end
    end
end
M_m=M_m.*K;
M_m(logical(eye(size(M_m))))=(1-d)*ones(size(M_m,1),1);
stability_m(1,k)=-log(max(abs(eig(M_m))));
k=k+1;
end


% 画图


subplot(1,4,4)
set(gca,'unit', 'centimeters', 'position', [20, 1, 5, 4]);

set(gca,'FontSize',8)
hold on 
stability=[mean(stability_r) mean(stability_e) mean(stability_c) mean(stability_m)];

b=bar(stability,0.4)
b.FaceColor = 'flat'; % 将 FaceColor 属性设置为 'flat' 以启用逐柱子颜色设置
b.CData(1, :) = [0 0.4470 0.7410]; % 第一个柱子的颜色
b.CData(2, :) = [0.8500 0.3250 0.0980]; % 第二个柱子的颜色
b.CData(3, :) = [0.9290 0.6940 0.1250]; % 第三个柱子的颜色
b.CData(4, :) = [0.4940 0.1840 0.5560]; % 第四个柱子的颜色

% set(b(1),'FaceColor',[0 0.4470 0.7410]);
% set(b(2),'FaceColor',[0.8500 0.3250 0.0980]);
% set(b(3),'FaceColor',[0.9290 0.6940 0.1250]);
% set(b(4),'FaceColor',[0.4940 0.1840 0.5560]);
set(b,'edgecolor','none');
%alpha(0.85);
axis([0.5 4.5 0 0.95]);
set(gca, 'xticklabel',{'Random','+/-','-/-','+/+'},'FontSize',8);
box on

ylabel('Stability','FontSize',9); 
xlabel('Interaction type','fontsize',9);
%title('Empirical dataset','fontsize',9);
title('$$d=1.5$$','Interpreter','latex','fontsize',9)