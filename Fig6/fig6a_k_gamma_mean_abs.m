

% This file can be used to generate Fig. 6a



%% 横坐标为k 纵坐标为E(|X|) 表示相互作用强度

%
clear all
clc



figure
set(gcf,'unit','centimeters','position',[15 5 10 6]) %长9 宽6
subplot(1,1,1)
set(gca,'unit', 'centimeters', 'position', [2 1 5 4]);
set(gca,'FontSize',8) %坐标轴字体的大小

%参数

k=0:0.01:6;
%E_abs=sigma*sqrt(k/(k+1));
E_abs=sqrt(k./(k+1));
plot(k,E_abs,'linestyle','-',"Color",'[0 0.4470 0.7410]','LineWidth',1)
hold on

axis ([0,5.3,0,1.05])


xlabel('$$k$$','Interpreter','latex','fontsize',9);
ylabel('The expected magnitude of the interactions','fontsize',9); 
%ylabel('$$E(|z|)/\sigma$$','Interpreter','latex','fontsize',9); 
box on