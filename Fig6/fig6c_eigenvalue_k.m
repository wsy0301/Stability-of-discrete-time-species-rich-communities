
% This file can be used to generate Fig. 6c 



%% =========================================================
% Random versus Exploitative communities
%
% Row 1:
%       d = 0.4
%
%       Random:
%           k = 0.5, 1, 2
%
%       Exploitative:
%           k = 0.5, 1, 2
%
%
% Row 2:
%       d = 1
%
%       Random:
%           k = 0.5, 1, 2
%
%       Exploitative:
%           k = 0.5, 1, 2
%
%
% 每个参数组合独立生成两个矩阵 realization
%
% 所有随机数都在构造矩阵过程中直接生成
%
% 不固定 topology / sign / interaction strength
%% =========================================================

clear all
clc
close all


%% =========================================================
% Parameters
%% =========================================================

sigma = 0.05;

S = 100;

C = 0.2;

d1=0.4;

d2=1;

k1=0.5;
k2=1;
k3=2;
%% =========================================================
% beta
%
% beta = sigma*sqrt(S*C)
%% =========================================================

beta = sigma*sqrt(S*C);


%% =========================================================
% Color
%% =========================================================

redColor = [228 26 28]/255;

%grayColor = [0.55 0.55 0.55];
grayColor=[31,120,180]/255;


equalColor = [0.30 0.30 0.30];


%% =========================================================
% Figure
%% =========================================================

figure

set(gcf,...
    'unit','centimeters',...
    'position',[2 2 26 10]);


%% =========================================================
% 公共变量
%% =========================================================

theta = linspace(0,2*pi,1500);


%% =========================================================
% Stable region |lambda| < 1
%% =========================================================

r = (0:0.01:1)';

theta_circle = pi*(-1:0.01:1);

x_circle = r*cos(theta_circle);

y_circle = r*sin(theta_circle);

circle_value = ...
    x_circle.^2 + y_circle.^2;


%% =========================================================
% Colormap
%% =========================================================

start_color = [253,174,107]/255;

end_color = [255,245,235]/255;

n_colors = 256;

custom_colormap = ...
    [linspace(start_color(1),end_color(1),n_colors)',...
     linspace(start_color(2),end_color(2),n_colors)',...
     linspace(start_color(3),end_color(3),n_colors)'];

colormap(custom_colormap);



%% =========================================================
%% PANEL 1
%
% Random
%
% d = 0.4
% k = 0.5
%% =========================================================

subplot(2,6,1)

set(gca,...
    'unit','centimeters',...
    'position',[1.0 5.5 3.5 3.5]);

hold on


d = d1;

k = k1;

c = 1-d;


shape = k;

scale = sigma*sqrt(1/(k^2+k));


% ---------------------------------------------------------
% Background unstable region
% ---------------------------------------------------------

patch(...
    [-1.1 1.1 1.1 -1.1],...
    [-1.1 -1.1 1.1 1.1],...
    [0.75 0.75 0.75],...
    'EdgeColor','none',...
    'FaceAlpha',0.5);


% ---------------------------------------------------------
% Stable region
% ---------------------------------------------------------

p_circle = pcolor(...
    x_circle,...
    y_circle,...
    4.5*(1-circle_value));

set(p_circle,'LineStyle','none');


% ---------------------------------------------------------
% Unit circle
% ---------------------------------------------------------

plot(...
    cos(theta),...
    sin(theta),...
    '-',...
    'Color',[0.35 0.35 0.35],...
    'LineWidth',0.8);


% ---------------------------------------------------------
% Numerical realization 1
% ---------------------------------------------------------

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;
                p3 = rand;


                if p2 <= 0.5

                    M(i,j) = gamrnd(shape,scale);

                else

                    M(i,j) = -gamrnd(shape,scale);

                end


                if p3 <= 0.5

                    M(j,i) = gamrnd(shape,scale);

                else

                    M(j,i) = -gamrnd(shape,scale);

                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% ---------------------------------------------------------
% Numerical realization 2
% ---------------------------------------------------------

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;
                p3 = rand;


                if p2 <= 0.5

                    M(i,j) = gamrnd(shape,scale);

                else

                    M(i,j) = -gamrnd(shape,scale);

                end


                if p3 <= 0.5

                    M(j,i) = gamrnd(shape,scale);

                else

                    M(j,i) = -gamrnd(shape,scale);

                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% ---------------------------------------------------------
% Random theoretical circle
% ---------------------------------------------------------

x_theory = c + beta*cos(theta);

y_theory = beta*sin(theta);


plot(...
    x_theory,...
    y_theory,...
    'k-',...
    'LineWidth',1);


% ---------------------------------------------------------
% Dominant eigenvalue
% ---------------------------------------------------------

lambda_dom = c + beta;

rho = abs(lambda_dom);


% equal-modulus circle
% plot(...
%     rho*cos(theta),...
%     rho*sin(theta),...
%     '--',...
%     'Color',equalColor,...
%     'LineWidth',0.8);


% dominant eigenvalue
plot(...
    real(lambda_dom),...
    imag(lambda_dom),...
    'o',...
    'MarkerSize',3,...
    'MarkerFaceColor',redColor,...
    'MarkerEdgeColor',redColor);


% ---------------------------------------------------------
% Axes
% ---------------------------------------------------------

% plot([-1.1 1.1],[0 0],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot([0 0],[-1.1 1.1],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot(0,0,'k+','MarkerSize',5);


axis equal

axis([-1.1 1.1 -1.1 1.1])

xticks([-1 0 1])

yticks([-1 0 1])

set(gca,...
    'FontSize',8,...
    'Layer','top');

box on


ylabel(...
    '$\mathrm{Im}(\lambda)$',...
    'Interpreter','latex',...
    'FontSize',9);


title(...
    {'Random','$k=0.5$'},...
    'Interpreter','latex',...
    'FontSize',8);



%% =========================================================
% PANEL 2
%
% Random
%
% d = 0.4
% k = 1
% =========================================================

subplot(2,6,2)

set(gca,...
    'unit','centimeters',...
    'position',[5.1 5.5 3.5 3.5]);

hold on


d = d1;

k = k2;

c = 1-d;


shape = k;

scale = sigma*sqrt(1/(k^2+k));


patch(...
    [-1.1 1.1 1.1 -1.1],...
    [-1.1 -1.1 1.1 1.1],...
    [0.75 0.75 0.75],...
    'EdgeColor','none',...
    'FaceAlpha',0.5);


p_circle = pcolor(...
    x_circle,...
    y_circle,...
    4.5*(1-circle_value));

set(p_circle,'LineStyle','none');


plot(...
    cos(theta),...
    sin(theta),...
    '-',...
    'Color',[0.35 0.35 0.35],...
    'LineWidth',0.8);


% realization 1

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;
                p3 = rand;


                if p2 <= 0.5
                    M(i,j) = gamrnd(shape,scale);
                else
                    M(i,j) = -gamrnd(shape,scale);
                end


                if p3 <= 0.5
                    M(j,i) = gamrnd(shape,scale);
                else
                    M(j,i) = -gamrnd(shape,scale);
                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% realization 2

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;
                p3 = rand;


                if p2 <= 0.5
                    M(i,j) = gamrnd(shape,scale);
                else
                    M(i,j) = -gamrnd(shape,scale);
                end


                if p3 <= 0.5
                    M(j,i) = gamrnd(shape,scale);
                else
                    M(j,i) = -gamrnd(shape,scale);
                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% theory

x_theory = c + beta*cos(theta);

y_theory = beta*sin(theta);


plot(...
    x_theory,...
    y_theory,...
    'k-',...
    'LineWidth',1);


lambda_dom = c + beta;

rho = abs(lambda_dom);


% plot(...
%     rho*cos(theta),...
%     rho*sin(theta),...
%     '--',...
%     'Color',equalColor,...
%     'LineWidth',0.8);


plot(...
    lambda_dom,...
    0,...
    'o',...
    'MarkerSize',3,...
    'MarkerFaceColor',redColor,...
    'MarkerEdgeColor',redColor);

% 
% plot([-1.1 1.1],[0 0],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot([0 0],[-1.1 1.1],...
%     '--','Color',[0.65 0.65 0.65]);

% plot(0,0,'k+','MarkerSize',5);


axis equal

axis([-1.1 1.1 -1.1 1.1])

xticks([-1 0 1])

yticks([-1 0 1])

set(gca,...
    'FontSize',8,...
    'Layer','top');

box on


title(...
    {'Random','$k=1$'},...
    'Interpreter','latex',...
    'FontSize',8);



%% =========================================================
%% PANEL 3
%
% Random
%
% d = 0.4
% k = 2
%% =========================================================

subplot(2,6,3)

set(gca,...
    'unit','centimeters',...
    'position',[9.2 5.5 3.5 3.5]);

hold on


d = d1;

k = k3;

c = 1-d;


shape = k;

scale = sigma*sqrt(1/(k^2+k));


patch(...
    [-1.1 1.1 1.1 -1.1],...
    [-1.1 -1.1 1.1 1.1],...
    [0.75 0.75 0.75],...
    'EdgeColor','none',...
    'FaceAlpha',0.5);


p_circle = pcolor(...
    x_circle,...
    y_circle,...
    4.5*(1-circle_value));

set(p_circle,'LineStyle','none');


plot(...
    cos(theta),...
    sin(theta),...
    '-',...
    'Color',[0.35 0.35 0.35],...
    'LineWidth',0.8);


% realization 1

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;
                p3 = rand;


                if p2 <= 0.5
                    M(i,j) = gamrnd(shape,scale);
                else
                    M(i,j) = -gamrnd(shape,scale);
                end


                if p3 <= 0.5
                    M(j,i) = gamrnd(shape,scale);
                else
                    M(j,i) = -gamrnd(shape,scale);
                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% realization 2

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;
                p3 = rand;


                if p2 <= 0.5
                    M(i,j) = gamrnd(shape,scale);
                else
                    M(i,j) = -gamrnd(shape,scale);
                end


                if p3 <= 0.5
                    M(j,i) = gamrnd(shape,scale);
                else
                    M(j,i) = -gamrnd(shape,scale);
                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% theory

x_theory = c + beta*cos(theta);

y_theory = beta*sin(theta);


plot(...
    x_theory,...
    y_theory,...
    'k-',...
    'LineWidth',1);


lambda_dom = c + beta;

rho = abs(lambda_dom);


% plot(...
%     rho*cos(theta),...
%     rho*sin(theta),...
%     '--',...
%     'Color',equalColor,...
%     'LineWidth',0.8);


plot(...
    lambda_dom,...
    0,...
    'o',...
    'MarkerSize',3,...
    'MarkerFaceColor',redColor,...
    'MarkerEdgeColor',redColor);


% plot([-1.1 1.1],[0 0],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot([0 0],[-1.1 1.1],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot(0,0,'k+','MarkerSize',5);


axis equal

axis([-1.1 1.1 -1.1 1.1])

xticks([-1 0 1])

yticks([-1 0 1])

set(gca,...
    'FontSize',8,...
    'Layer','top');

box on


title(...
    {'Random','$k=2$'},...
    'Interpreter','latex',...
    'FontSize',8);



%% =========================================================
% PANEL 4
%
% Exploitative
%
% d = 0.4
% k = 0.5
% =========================================================

subplot(2,6,4)

set(gca,...
    'unit','centimeters',...
    'position',[13.6 5.5 3.5 3.5]);

hold on


d = d1;

k = k1;

c = 1-d;


shape = k;

scale = sigma*sqrt(1/(k^2+k));


patch(...
    [-1.1 1.1 1.1 -1.1],...
    [-1.1 -1.1 1.1 1.1],...
    [0.75 0.75 0.75],...
    'EdgeColor','none',...
    'FaceAlpha',0.5);


p_circle = pcolor(...
    x_circle,...
    y_circle,...
    4.5*(1-circle_value));

set(p_circle,'LineStyle','none');


plot(...
    cos(theta),...
    sin(theta),...
    '-',...
    'Color',[0.35 0.35 0.35],...
    'LineWidth',0.8);


% realization 1

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;


                if p2 <= 0.5

                    M(i,j) = gamrnd(shape,scale);

                    M(j,i) = -gamrnd(shape,scale);

                else

                    M(i,j) = -gamrnd(shape,scale);

                    M(j,i) = gamrnd(shape,scale);

                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% realization 2

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;


                if p2 <= 0.5

                    M(i,j) = gamrnd(shape,scale);

                    M(j,i) = -gamrnd(shape,scale);

                else

                    M(i,j) = -gamrnd(shape,scale);

                    M(j,i) = gamrnd(shape,scale);

                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% theoretical ellipse

a = beta/(k+1);

b = beta*(2*k+1)/(k+1);


x_ell = c + a*cos(theta);

y_ell = b*sin(theta);


plot(...
    x_ell,...
    y_ell,...
    'k-',...
    'LineWidth',1);


% dominant eigenvalue

lambda_right = c+a;

lambda_left = c-a;


rho_right = abs(lambda_right);

rho_left = abs(lambda_left);


if b > a

    u_dom = a*c/(b^2-a^2);

else

    u_dom = 2;

end


if abs(u_dom) <= 1

    theta_dom = acos(u_dom);

    lambda_arc_up = ...
        c ...
        + a*cos(theta_dom) ...
        + 1i*b*sin(theta_dom);

    lambda_arc_down = conj(lambda_arc_up);

    rho_arc = abs(lambda_arc_up);

else

    rho_arc = -Inf;

end


rho_endpoint = max(rho_right,rho_left);


if rho_arc > rho_endpoint

    rho = rho_arc;

% 
%     plot(...
%         rho*cos(theta),...
%         rho*sin(theta),...
%         '--',...
%         'Color',equalColor,...
%         'LineWidth',0.8);


    plot(...
        real(lambda_arc_up),...
        imag(lambda_arc_up),...
        'o',...
        'MarkerSize',3,...
        'MarkerFaceColor',redColor,...
        'MarkerEdgeColor',redColor);


    plot(...
        real(lambda_arc_down),...
        imag(lambda_arc_down),...
        'o',...
        'MarkerSize',3,...
        'MarkerFaceColor',redColor,...
        'MarkerEdgeColor',redColor);

else

    if rho_right >= rho_left

        lambda_dom = lambda_right;

    else

        lambda_dom = lambda_left;

    end


    rho = abs(lambda_dom);

% 
%     plot(...
%         rho*cos(theta),...
%         rho*sin(theta),...
%         '--',...
%         'Color',equalColor,...
%         'LineWidth',0.8);


    plot(...
        real(lambda_dom),...
        imag(lambda_dom),...
        'o',...
        'MarkerSize',3,...
        'MarkerFaceColor',redColor,...
        'MarkerEdgeColor',redColor);

end

% 
% plot([-1.1 1.1],[0 0],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot([0 0],[-1.1 1.1],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot(0,0,'k+','MarkerSize',5);


axis equal

axis([-1.1 1.1 -1.1 1.1])

xticks([-1 0 1])

yticks([-1 0 1])

set(gca,...
    'FontSize',8,...
    'Layer','top');

box on


title(...
    {'Exploitative','$k=0.5$'},...
    'Interpreter','latex',...
    'FontSize',8);



%% =========================================================
% PANEL 5
%
% Exploitative
%
% d = 0.4
% k = 1
% =========================================================

subplot(2,6,5)

set(gca,...
    'unit','centimeters',...
    'position',[17.7 5.5 3.5 3.5]);

hold on


d = d1;

k = k2;

c = 1-d;


shape = k;

scale = sigma*sqrt(1/(k^2+k));


patch(...
    [-1.1 1.1 1.1 -1.1],...
    [-1.1 -1.1 1.1 1.1],...
    [0.75 0.75 0.75],...
    'EdgeColor','none',...
    'FaceAlpha',0.5);


p_circle = pcolor(...
    x_circle,...
    y_circle,...
    4.5*(1-circle_value));

set(p_circle,'LineStyle','none');


plot(...
    cos(theta),...
    sin(theta),...
    '-',...
    'Color',[0.35 0.35 0.35],...
    'LineWidth',0.8);


% realization 1

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;


                if p2 <= 0.5

                    M(i,j) = gamrnd(shape,scale);

                    M(j,i) = -gamrnd(shape,scale);

                else

                    M(i,j) = -gamrnd(shape,scale);

                    M(j,i) = gamrnd(shape,scale);

                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% realization 2

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;


                if p2 <= 0.5

                    M(i,j) = gamrnd(shape,scale);

                    M(j,i) = -gamrnd(shape,scale);

                else

                    M(i,j) = -gamrnd(shape,scale);

                    M(j,i) = gamrnd(shape,scale);

                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% ellipse

a = beta/(k+1);

b = beta*(2*k+1)/(k+1);


x_ell = c+a*cos(theta);

y_ell = b*sin(theta);


plot(...
    x_ell,...
    y_ell,...
    'k-',...
    'LineWidth',1);


% dominant

lambda_right = c+a;

lambda_left = c-a;

rho_right = abs(lambda_right);

rho_left = abs(lambda_left);


u_dom = a*c/(b^2-a^2);


if abs(u_dom) <= 1

    theta_dom = acos(u_dom);

    lambda_arc_up = ...
        c+a*cos(theta_dom) ...
        +1i*b*sin(theta_dom);

    lambda_arc_down = conj(lambda_arc_up);

    rho_arc = abs(lambda_arc_up);

else

    rho_arc = -Inf;

end


rho_endpoint = max(rho_right,rho_left);


if rho_arc > rho_endpoint

    rho = rho_arc;


%     plot(...
%         rho*cos(theta),...
%         rho*sin(theta),...
%         '--',...
%         'Color',equalColor,...
%         'LineWidth',0.8);


    plot(...
        real(lambda_arc_up),...
        imag(lambda_arc_up),...
        'o',...
        'MarkerSize',3,...
        'MarkerFaceColor',redColor,...
        'MarkerEdgeColor',redColor);


    plot(...
        real(lambda_arc_down),...
        imag(lambda_arc_down),...
        'o',...
        'MarkerSize',3,...
        'MarkerFaceColor',redColor,...
        'MarkerEdgeColor',redColor);

else

    if rho_right >= rho_left
        lambda_dom = lambda_right;
    else
        lambda_dom = lambda_left;
    end


    rho = abs(lambda_dom);


%     plot(...
%         rho*cos(theta),...
%         rho*sin(theta),...
%         '--',...
%         'Color',equalColor,...
%         'LineWidth',0.8);


    plot(...
        real(lambda_dom),...
        imag(lambda_dom),...
        'o',...
        'MarkerSize',3,...
        'MarkerFaceColor',redColor,...
        'MarkerEdgeColor',redColor);

end


% plot([-1.1 1.1],[0 0],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot([0 0],[-1.1 1.1],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot(0,0,'k+','MarkerSize',5);


axis equal

axis([-1.1 1.1 -1.1 1.1])

xticks([-1 0 1])

yticks([-1 0 1])

set(gca,...
    'FontSize',8,...
    'Layer','top');

box on


title(...
    {'Exploitative','$k=1$'},...
    'Interpreter','latex',...
    'FontSize',8);



%% =========================================================
% PANEL 6
%
% Exploitative
%
% d = 0.4
% k = 2
%% =========================================================

subplot(2,6,6)

set(gca,...
    'unit','centimeters',...
    'position',[21.8 5.5 3.5 3.5]);

hold on


d = d1;

k = k3;

c = 1-d;


shape = k;

scale = sigma*sqrt(1/(k^2+k));


patch(...
    [-1.1 1.1 1.1 -1.1],...
    [-1.1 -1.1 1.1 1.1],...
    [0.75 0.75 0.75],...
    'EdgeColor','none',...
    'FaceAlpha',0.5);


p_circle = pcolor(...
    x_circle,...
    y_circle,...
    4.5*(1-circle_value));

set(p_circle,'LineStyle','none');


plot(...
    cos(theta),...
    sin(theta),...
    '-',...
    'Color',[0.35 0.35 0.35],...
    'LineWidth',0.8);


% realization 1

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;


                if p2 <= 0.5

                    M(i,j) = gamrnd(shape,scale);

                    M(j,i) = -gamrnd(shape,scale);

                else

                    M(i,j) = -gamrnd(shape,scale);

                    M(j,i) = gamrnd(shape,scale);

                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% realization 2

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;


                if p2 <= 0.5

                    M(i,j) = gamrnd(shape,scale);

                    M(j,i) = -gamrnd(shape,scale);

                else

                    M(i,j) = -gamrnd(shape,scale);

                    M(j,i) = gamrnd(shape,scale);

                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% ellipse

a = beta/(k+1);

b = beta*(2*k+1)/(k+1);


x_ell = c+a*cos(theta);

y_ell = b*sin(theta);


plot(...
    x_ell,...
    y_ell,...
    'k-',...
    'LineWidth',1);


% dominant

lambda_right = c+a;

lambda_left = c-a;

rho_right = abs(lambda_right);

rho_left = abs(lambda_left);


u_dom = a*c/(b^2-a^2);


if abs(u_dom) <= 1

    theta_dom = acos(u_dom);

    lambda_arc_up = ...
        c+a*cos(theta_dom) ...
        +1i*b*sin(theta_dom);

    lambda_arc_down = conj(lambda_arc_up);

    rho_arc = abs(lambda_arc_up);

else

    rho_arc = -Inf;

end


rho_endpoint = max(rho_right,rho_left);


if rho_arc > rho_endpoint

    rho = rho_arc;


%     plot(...
%         rho*cos(theta),...
%         rho*sin(theta),...
%         '--',...
%         'Color',equalColor,...
%         'LineWidth',0.8);


    plot(...
        real(lambda_arc_up),...
        imag(lambda_arc_up),...
        'o',...
        'MarkerSize',3,...
        'MarkerFaceColor',redColor,...
        'MarkerEdgeColor',redColor);


    plot(...
        real(lambda_arc_down),...
        imag(lambda_arc_down),...
        'o',...
        'MarkerSize',3,...
        'MarkerFaceColor',redColor,...
        'MarkerEdgeColor',redColor);

else

    if rho_right >= rho_left
        lambda_dom = lambda_right;
    else
        lambda_dom = lambda_left;
    end


    rho = abs(lambda_dom);


%     plot(...
%         rho*cos(theta),...
%         rho*sin(theta),...
%         '--',...
%         'Color',equalColor,...
%         'LineWidth',0.8);


    plot(...
        real(lambda_dom),...
        imag(lambda_dom),...
        'o',...
        'MarkerSize',5,...
        'MarkerFaceColor',redColor,...
        'MarkerEdgeColor',redColor);

end


% plot([-1.1 1.1],[0 0],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot([0 0],[-1.1 1.1],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot(0,0,'k+','MarkerSize',5);


axis equal

axis([-1.1 1.1 -1.1 1.1])

xticks([-1 0 1])

yticks([-1 0 1])

set(gca,...
    'FontSize',8,...
    'Layer','top');

box on


title(...
    {'Exploitative','$k=2$'},...
    'Interpreter','latex',...
    'FontSize',8);



%% =========================================================
% PANEL 7
%
% Random
%
% d = 1
% k = 0.5
%% =========================================================

subplot(2,6,7)

set(gca,...
    'unit','centimeters',...
    'position',[1.0 1.0 3.5 3.5]);

hold on


d = d2;

k = k1;

c = 1-d;


shape = k;

scale = sigma*sqrt(1/(k^2+k));


patch(...
    [-1.1 1.1 1.1 -1.1],...
    [-1.1 -1.1 1.1 1.1],...
    [0.75 0.75 0.75],...
    'EdgeColor','none',...
    'FaceAlpha',0.5);


p_circle = pcolor(...
    x_circle,...
    y_circle,...
    4.5*(1-circle_value));

set(p_circle,'LineStyle','none');


plot(...
    cos(theta),...
    sin(theta),...
    '-',...
    'Color',[0.35 0.35 0.35],...
    'LineWidth',0.8);


% realization 1

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;
                p3 = rand;


                if p2 <= 0.5
                    M(i,j) = gamrnd(shape,scale);
                else
                    M(i,j) = -gamrnd(shape,scale);
                end


                if p3 <= 0.5
                    M(j,i) = gamrnd(shape,scale);
                else
                    M(j,i) = -gamrnd(shape,scale);
                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% realization 2

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;
                p3 = rand;


                if p2 <= 0.5
                    M(i,j) = gamrnd(shape,scale);
                else
                    M(i,j) = -gamrnd(shape,scale);
                end


                if p3 <= 0.5
                    M(j,i) = gamrnd(shape,scale);
                else
                    M(j,i) = -gamrnd(shape,scale);
                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% theoretical circle
%
% d=1:
%
% center = 0
%
% 圆周上的理论边界点模长全部相同
% ---------------------------------------------------------

% theoretical random circle
% theory

x_theory = c + beta*cos(theta);

y_theory = beta*sin(theta);


plot(...
    x_theory,...
    y_theory,...
    'k-',...
    'LineWidth',1);


lambda_dom = c + beta;

% plot(...
%     rho*cos(theta),...
%     rho*sin(theta),...
%     '--',...
%     'Color',equalColor,...
%     'LineWidth',0.8);


plot(...
    lambda_dom,...
    0,...
    'o',...
    'MarkerSize',3,...
    'MarkerFaceColor',redColor,...
    'MarkerEdgeColor',redColor);

% plot([-1.1 1.1],[0 0],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot([0 0],[-1.1 1.1],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot(0,0,'k+','MarkerSize',5);


axis equal

axis([-1.1 1.1 -1.1 1.1])

xticks([-1 0 1])

yticks([-1 0 1])

set(gca,...
    'FontSize',8,...
    'Layer','top');

box on


xlabel(...
    '$\mathrm{Re}(\lambda)$',...
    'Interpreter','latex',...
    'FontSize',9);

ylabel(...
    '$\mathrm{Im}(\lambda)$',...
    'Interpreter','latex',...
    'FontSize',9);


title(...
    {'Random','$k=0.5$'},...
    'Interpreter','latex',...
    'FontSize',8);



%% =========================================================
% PANEL 8
%
% Random
%
% d = 1
% k = 1
% =========================================================

subplot(2,6,8)

set(gca,...
    'unit','centimeters',...
    'position',[5.1 1.0 3.5 3.5]);

hold on


d = d2;

k = k2;

c = 1-d;


shape = k;

scale = sigma*sqrt(1/(k^2+k));


patch(...
    [-1.1 1.1 1.1 -1.1],...
    [-1.1 -1.1 1.1 1.1],...
    [0.75 0.75 0.75],...
    'EdgeColor','none',...
    'FaceAlpha',0.5);


p_circle = pcolor(...
    x_circle,...
    y_circle,...
    4.5*(1-circle_value));

set(p_circle,'LineStyle','none');


plot(...
    cos(theta),...
    sin(theta),...
    '-',...
    'Color',[0.35 0.35 0.35],...
    'LineWidth',0.8);


% realization 1

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;
                p3 = rand;


                if p2 <= 0.5
                    M(i,j) = gamrnd(shape,scale);
                else
                    M(i,j) = -gamrnd(shape,scale);
                end


                if p3 <= 0.5
                    M(j,i) = gamrnd(shape,scale);
                else
                    M(j,i) = -gamrnd(shape,scale);
                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% realization 2

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;
                p3 = rand;


                if p2 <= 0.5
                    M(i,j) = gamrnd(shape,scale);
                else
                    M(i,j) = -gamrnd(shape,scale);
                end


                if p3 <= 0.5
                    M(j,i) = gamrnd(shape,scale);
                else
                    M(j,i) = -gamrnd(shape,scale);
                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% theoretical circle
%
% d=1:
%
% center = 0
%
% 圆周上的理论边界点模长全部相同
% ---------------------------------------------------------

% theoretical random circle
% theory

x_theory = c + beta*cos(theta);

y_theory = beta*sin(theta);


plot(...
    x_theory,...
    y_theory,...
    'k-',...
    'LineWidth',1);


lambda_dom = c + beta;

% plot(...
%     rho*cos(theta),...
%     rho*sin(theta),...
%     '--',...
%     'Color',equalColor,...
%     'LineWidth',0.8);


plot(...
    lambda_dom,...
    0,...
    'o',...
    'MarkerSize',3,...
    'MarkerFaceColor',redColor,...
    'MarkerEdgeColor',redColor);


axis equal

axis([-1.1 1.1 -1.1 1.1])

xticks([-1 0 1])

yticks([-1 0 1])

set(gca,...
    'FontSize',8,...
    'Layer','top');

box on


xlabel(...
    '$\mathrm{Re}(\lambda)$',...
    'Interpreter','latex',...
    'FontSize',9);


title(...
    {'Random','$k=1$'},...
    'Interpreter','latex',...
    'FontSize',8);



%% =========================================================
%% PANEL 9
%
% Random
%
% d = 1
% k = 2
%% =========================================================

subplot(2,6,9)

set(gca,...
    'unit','centimeters',...
    'position',[9.2 1.0 3.5 3.5]);

hold on


d = d2;

k = k3;

c = 1-d;


shape = k;

scale = sigma*sqrt(1/(k^2+k));


patch(...
    [-1.1 1.1 1.1 -1.1],...
    [-1.1 -1.1 1.1 1.1],...
    [0.75 0.75 0.75],...
    'EdgeColor','none',...
    'FaceAlpha',0.5);


p_circle = pcolor(...
    x_circle,...
    y_circle,...
    4.5*(1-circle_value));

set(p_circle,'LineStyle','none');


plot(...
    cos(theta),...
    sin(theta),...
    '-',...
    'Color',[0.35 0.35 0.35],...
    'LineWidth',0.8);


% realization 1

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;
                p3 = rand;


                if p2 <= 0.5
                    M(i,j) = gamrnd(shape,scale);
                else
                    M(i,j) = -gamrnd(shape,scale);
                end


                if p3 <= 0.5
                    M(j,i) = gamrnd(shape,scale);
                else
                    M(j,i) = -gamrnd(shape,scale);
                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% realization 2

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 1-d;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;
                p3 = rand;


                if p2 <= 0.5
                    M(i,j) = gamrnd(shape,scale);
                else
                    M(i,j) = -gamrnd(shape,scale);
                end


                if p3 <= 0.5
                    M(j,i) = gamrnd(shape,scale);
                else
                    M(j,i) = -gamrnd(shape,scale);
                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% theoretical circle
%
% d=1:
%
% center = 0
%
% 圆周上的理论边界点模长全部相同
% ---------------------------------------------------------

% theoretical random circle
% theory

x_theory = c + beta*cos(theta);

y_theory = beta*sin(theta);


plot(...
    x_theory,...
    y_theory,...
    'k-',...
    'LineWidth',1);


lambda_dom = c + beta;

% plot(...
%     rho*cos(theta),...
%     rho*sin(theta),...
%     '--',...
%     'Color',equalColor,...
%     'LineWidth',0.8);


plot(...
    lambda_dom,...
    0,...
    'o',...
    'MarkerSize',3,...
    'MarkerFaceColor',redColor,...
    'MarkerEdgeColor',redColor);

axis equal

axis([-1.1 1.1 -1.1 1.1])

xticks([-1 0 1])

yticks([-1 0 1])

set(gca,...
    'FontSize',8,...
    'Layer','top');

box on


xlabel(...
    '$\mathrm{Re}(\lambda)$',...
    'Interpreter','latex',...
    'FontSize',9);


title(...
    {'Random','$k=2$'},...
    'Interpreter','latex',...
    'FontSize',8);



%% =========================================================
%% PANEL 10
%
% Exploitative
%
% d = 1
% k = 0.5
%% =========================================================

subplot(2,6,10)

set(gca,...
    'unit','centimeters',...
    'position',[13.6 1.0 3.5 3.5]);

hold on


d = d2;

k = k1;

c = 1-d;


shape = k;

scale = sigma*sqrt(1/(k^2+k));


patch(...
    [-1.1 1.1 1.1 -1.1],...
    [-1.1 -1.1 1.1 1.1],...
    [0.75 0.75 0.75],...
    'EdgeColor','none',...
    'FaceAlpha',0.5);


p_circle = pcolor(...
    x_circle,...
    y_circle,...
    4.5*(1-circle_value));

set(p_circle,'LineStyle','none');


plot(...
    cos(theta),...
    sin(theta),...
    '-',...
    'Color',[0.35 0.35 0.35],...
    'LineWidth',0.8);


% realization 1

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 0;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;


                if p2 <= 0.5

                    M(i,j) = gamrnd(shape,scale);

                    M(j,i) = -gamrnd(shape,scale);

                else

                    M(i,j) = -gamrnd(shape,scale);

                    M(j,i) = gamrnd(shape,scale);

                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% realization 2

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 0;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;


                if p2 <= 0.5

                    M(i,j) = gamrnd(shape,scale);

                    M(j,i) = -gamrnd(shape,scale);

                else

                    M(i,j) = -gamrnd(shape,scale);

                    M(j,i) = gamrnd(shape,scale);

                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% ellipse

a = beta/(k+1);

b = beta*(2*k+1)/(k+1);


plot(...
    a*cos(theta),...
    b*sin(theta),...
    'k-',...
    'LineWidth',1);


% d = 1
%
% b > a
%
% vertical endpoints are dominant
%% ---------------------------------------------------------

% rho = b;
% 
% 
% plot(...
%     rho*cos(theta),...
%     rho*sin(theta),...
%     '--',...
%     'Color',equalColor,...
%     'LineWidth',0.8);


plot(...
    0,b,...
    'o',...
    'MarkerSize',3,...
    'MarkerFaceColor',redColor,...
    'MarkerEdgeColor',redColor);


plot(...
    0,-b,...
    'o',...
    'MarkerSize',3,...
    'MarkerFaceColor',redColor,...
    'MarkerEdgeColor',redColor);


% plot([-1.1 1.1],[0 0],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot([0 0],[-1.1 1.1],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot(0,0,'k+','MarkerSize',5);


axis equal

axis([-1.1 1.1 -1.1 1.1])

xticks([-1 0 1])

yticks([-1 0 1])

set(gca,...
    'FontSize',8,...
    'Layer','top');

box on


xlabel(...
    '$\mathrm{Re}(\lambda)$',...
    'Interpreter','latex',...
    'FontSize',9);


title(...
    {'Exploitative','$k=0.5$'},...
    'Interpreter','latex',...
    'FontSize',8);



%% =========================================================
%% PANEL 11
%
% Exploitative
%
% d = 1
% k = 1
% =========================================================

subplot(2,6,11)

set(gca,...
    'unit','centimeters',...
    'position',[17.7 1.0 3.5 3.5]);

hold on


d = d2;

k = k2;


shape = k;

scale = sigma*sqrt(1/(k^2+k));


patch(...
    [-1.1 1.1 1.1 -1.1],...
    [-1.1 -1.1 1.1 1.1],...
    [0.75 0.75 0.75],...
    'EdgeColor','none',...
    'FaceAlpha',0.5);


p_circle = pcolor(...
    x_circle,...
    y_circle,...
    4.5*(1-circle_value));

set(p_circle,'LineStyle','none');


plot(...
    cos(theta),...
    sin(theta),...
    '-',...
    'Color',[0.35 0.35 0.35],...
    'LineWidth',0.8);


% realization 1

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 0;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;


                if p2 <= 0.5

                    M(i,j) = gamrnd(shape,scale);

                    M(j,i) = -gamrnd(shape,scale);

                else

                    M(i,j) = -gamrnd(shape,scale);

                    M(j,i) = gamrnd(shape,scale);

                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% realization 2

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 0;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;


                if p2 <= 0.5

                    M(i,j) = gamrnd(shape,scale);

                    M(j,i) = -gamrnd(shape,scale);

                else

                    M(i,j) = -gamrnd(shape,scale);

                    M(j,i) = gamrnd(shape,scale);

                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


a = beta/(k+1);

b = beta*(2*k+1)/(k+1);


plot(...
    a*cos(theta),...
    b*sin(theta),...
    'k-',...
    'LineWidth',1);


% rho = b;
% 
% 
% plot(...
%     rho*cos(theta),...
%     rho*sin(theta),...
%     '--',...
%     'Color',equalColor,...
%     'LineWidth',0.8);


plot(...
    0,b,...
    'o',...
    'MarkerSize',3,...
    'MarkerFaceColor',redColor,...
    'MarkerEdgeColor',redColor);


plot(...
    0,-b,...
    'o',...
    'MarkerSize',3,...
    'MarkerFaceColor',redColor,...
    'MarkerEdgeColor',redColor);


% plot([-1.1 1.1],[0 0],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot([0 0],[-1.1 1.1],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot(0,0,'k+','MarkerSize',5);


axis equal

axis([-1.1 1.1 -1.1 1.1])

xticks([-1 0 1])

yticks([-1 0 1])

set(gca,...
    'FontSize',8,...
    'Layer','top');

box on


xlabel(...
    '$\mathrm{Re}(\lambda)$',...
    'Interpreter','latex',...
    'FontSize',9);


title(...
    {'Exploitative','$k=1$'},...
    'Interpreter','latex',...
    'FontSize',8);



%% =========================================================
% PANEL 12
%
% Exploitative
%
% d = 1
% k = 2
% =========================================================

subplot(2,6,12)

set(gca,...
    'unit','centimeters',...
    'position',[21.8 1.0 3.5 3.5]);

hold on


d = 1;

k = k3;


shape = k;

scale = sigma*sqrt(1/(k^2+k));


patch(...
    [-1.1 1.1 1.1 -1.1],...
    [-1.1 -1.1 1.1 1.1],...
    [0.75 0.75 0.75],...
    'EdgeColor','none',...
    'FaceAlpha',0.5);


p_circle = pcolor(...
    x_circle,...
    y_circle,...
    4.5*(1-circle_value));

set(p_circle,'LineStyle','none');


plot(...
    cos(theta),...
    sin(theta),...
    '-',...
    'Color',[0.35 0.35 0.35],...
    'LineWidth',0.8);


% realization 1

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 0;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;


                if p2 <= 0.5

                    M(i,j) = gamrnd(shape,scale);

                    M(j,i) = -gamrnd(shape,scale);

                else

                    M(i,j) = -gamrnd(shape,scale);

                    M(j,i) = gamrnd(shape,scale);

                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


% realization 2

M = zeros(S);


for i = 1:S

    for j = 1:S

        if i == j

            M(i,j) = 0;

        elseif i > j

            p1 = rand;

            if p1 <= C

                p2 = rand;


                if p2 <= 0.5

                    M(i,j) = gamrnd(shape,scale);

                    M(j,i) = -gamrnd(shape,scale);

                else

                    M(i,j) = -gamrnd(shape,scale);

                    M(j,i) = gamrnd(shape,scale);

                end

            end

        end

    end

end


lambda_num = eig(M);


plot(...
    real(lambda_num),...
    imag(lambda_num),...
    '.',...
    'Color',grayColor,...
    'MarkerSize',5);


a = beta/(k+1);

b = beta*(2*k+1)/(k+1);


plot(...
    a*cos(theta),...
    b*sin(theta),...
    'k-',...
    'LineWidth',1);

% 
% rho = b;
% 
% 
% plot(...
%     rho*cos(theta),...
%     rho*sin(theta),...
%     '--',...
%     'Color',equalColor,...
%     'LineWidth',0.8);


plot(...
    0,b,...
    'o',...
    'MarkerSize',3,...
    'MarkerFaceColor',redColor,...
    'MarkerEdgeColor',redColor);


plot(...
    0,-b,...
    'o',...
    'MarkerSize',3,...
    'MarkerFaceColor',redColor,...
    'MarkerEdgeColor',redColor);


% plot([-1.1 1.1],[0 0],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot([0 0],[-1.1 1.1],...
%     '--','Color',[0.65 0.65 0.65]);
% 
% plot(0,0,'k+','MarkerSize',5);


axis equal

axis([-1.1 1.1 -1.1 1.1])

xticks([-1 0 1])

yticks([-1 0 1])

set(gca,...
    'FontSize',8,...
    'Layer','top');

box on


xlabel(...
    '$\mathrm{Re}(\lambda)$',...
    'Interpreter','latex',...
    'FontSize',9);


title(...
    {'Exploitative','$k=2$'},...
    'Interpreter','latex',...
    'FontSize',8);



%% =========================================================
% Output
%% =========================================================

fprintf('\n');

fprintf('beta = %.6f\n',beta);

fprintf('\n');


fprintf('Random theoretical radius:\n');

fprintf('beta = %.6f\n',beta);

fprintf('\n');


fprintf('Exploitative theoretical semi-axes:\n');


k = 0.5;

a = beta/(k+1);

b = beta*(2*k+1)/(k+1);

fprintf(...
    'k = 0.5: a = %.6f, b = %.6f\n',...
    a,b);


k = 1;

a = beta/(k+1);

b = beta*(2*k+1)/(k+1);

fprintf(...
    'k = 1:   a = %.6f, b = %.6f\n',...
    a,b);


k = 2;

a = beta/(k+1);

b = beta*(2*k+1)/(k+1);

fprintf(...
    'k = 2:   a = %.6f, b = %.6f\n',...
    a,b);