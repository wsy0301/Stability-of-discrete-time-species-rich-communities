
% This file can be used to generate Fig. 4d 



%% =========================================================
% Single-species Beverton-Holt model
%
% X(k+1) = r*X(k)/(1+s*X(k))
%
% Prescribed equilibrium:
%
%       X* = 1
%
% Therefore:
%
%       r = 1 + s
%
% Self-regulation strength:
%
%       d = s*X*/r
%
% Since X*=1:
%
%       d = s/(1+s)
%
%
% Community matrix:
%
%       M = r/(1+sX*)^2
%
% Since r=1+s and X*=1:
%
%       M = 1/(1+s)
%         = 1-d
%
%
% Stability:
%
%       rho = |M|
%
%       alpha = -log(rho)
%
%
% NEW recovery criterion:
%
% Initial perturbation:
%
%       D0 = |X(0)-X*|
%
% At time k:
%
%       D(k) = |X(k)-X*|
%
% Recovered when:
%
%       D(k) <= 0.005*D0
%
% i.e. the perturbation has declined to
% 0.5% of its initial magnitude.
%
%% =========================================================

clear all
clc
close all



%% =========================================================
% Common parameters
%% =========================================================

X_star = 1;


% Initial abundance
X0 = 0.8;


% Extinction threshold
extinction_threshold = 0.001;


% Maximum simulation time
max_step = 30;


k = 0:max_step;



%% =========================================================
% NEW recovery criterion
%% =========================================================

% Initial Euclidean distance from equilibrium
%
% For one species, Euclidean distance is simply:
%
%       |X0-X*|

D0 = ...
    abs(X0-X_star);


% 0.5% of initial perturbation
recovery_fraction = ...
    0.005;


% Recovery threshold
recovery_threshold = ...
    recovery_fraction*D0;



fprintf('\n')

disp('=============================================')
disp('Recovery criterion')
disp('=============================================')


fprintf(...
    'Initial perturbation D0 = %.6f\n',...
    D0);


fprintf(...
    'Recovery threshold = %.6f\n',...
    recovery_threshold);


fprintf(...
    'Threshold = %.2f%% of initial perturbation\n',...
    100*recovery_fraction);



%% =========================================================
%% Case 1
%
% s = 0.4
%% =========================================================

s1 = 0.4;


% Because X*=1:
r1 = 1 + s1;


% Feedback-matrix self-regulation
d1 = ...
    s1*X_star/r1;



%% =========================================================
% Community matrix
%% =========================================================

M1 = ...
    r1/(1+s1*X_star)^2;



%% =========================================================
% Dominant eigenvalue modulus
%% =========================================================

eig1 = ...
    eig(M1);


rho1 = ...
    max(abs(eig1));



%% =========================================================
% Stability
%% =========================================================

stab1 = ...
    -log(rho1);



%% =========================================================
% Simulation
%% =========================================================

X1 = ...
    zeros(1,max_step+1);


X1(1) = ...
    X0;



for ii = 1:max_step


    if X1(ii) == 0


        X1(ii+1) = 0;


    else


        X1(ii+1) = ...
            r1*X1(ii) ...
            /(1+s1*X1(ii));


        % Extinction criterion
        if X1(ii+1) < extinction_threshold


            X1(ii+1) = 0;


        end


    end


end



%% =========================================================
% Euclidean distance from equilibrium
%% =========================================================

D1 = ...
    abs(X1-X_star);



%% =========================================================
% Recovery time
%
% First k satisfying:
%
%       D(k) <= 0.005 D0
%% =========================================================

recovery1 = ...
    NaN;


for ii = 1:length(X1)


    if D1(ii) <= recovery_threshold


        recovery1 = ...
            k(ii);


        break


    end


end



%% =========================================================
%% Case 2
%
% s = 0.9
%% =========================================================

s2 = 0.9;


r2 = ...
    1 + s2;


d2 = ...
    s2*X_star/r2;



%% =========================================================
% Community matrix
%% =========================================================

M2 = ...
    r2/(1+s2*X_star)^2;



%% =========================================================
% Dominant eigenvalue
%% =========================================================

eig2 = ...
    eig(M2);


rho2 = ...
    max(abs(eig2));



%% =========================================================
% Stability
%% =========================================================

stab2 = ...
    -log(rho2);



%% =========================================================
% Simulation
%% =========================================================

X2 = ...
    zeros(1,max_step+1);


X2(1) = ...
    X0;



for ii = 1:max_step


    if X2(ii) == 0


        X2(ii+1) = 0;


    else


        X2(ii+1) = ...
            r2*X2(ii) ...
            /(1+s2*X2(ii));


        if X2(ii+1) < extinction_threshold


            X2(ii+1) = 0;


        end


    end


end



%% =========================================================
% Distance
%% =========================================================

D2 = ...
    abs(X2-X_star);



%% =========================================================
% Recovery time
%% =========================================================

recovery2 = ...
    NaN;


for ii = 1:length(X2)


    if D2(ii) <= recovery_threshold


        recovery2 = ...
            k(ii);


        break


    end


end



%% =========================================================
%% Case 3
%
% s = 1.73
%% =========================================================

s3 = 1.73;


r3 = ...
    1 + s3;


d3 = ...
    s3*X_star/r3;



%% =========================================================
% Community matrix
%% =========================================================

M3 = ...
    r3/(1+s3*X_star)^2;



%% =========================================================
% Dominant eigenvalue
%% =========================================================

eig3 = ...
    eig(M3);


rho3 = ...
    max(abs(eig3));



%% =========================================================
% Stability
%% =========================================================

stab3 = ...
    -log(rho3);



%% =========================================================
% Simulation
%% =========================================================

X3 = ...
    zeros(1,max_step+1);


X3(1) = ...
    X0;



for ii = 1:max_step


    if X3(ii) == 0


        X3(ii+1) = 0;


    else


        X3(ii+1) = ...
            r3*X3(ii) ...
            /(1+s3*X3(ii));


        if X3(ii+1) < extinction_threshold


            X3(ii+1) = 0;


        end


    end


end



%% =========================================================
% Distance
%% =========================================================

D3 = ...
    abs(X3-X_star);



%% =========================================================
% Recovery time
%% =========================================================

recovery3 = ...
    NaN;


for ii = 1:length(X3)


    if D3(ii) <= recovery_threshold


        recovery3 = ...
            k(ii);


        break


    end


end



%% =========================================================
%% Case 4
%
% s = 4
%% =========================================================

s4 = 4;


r4 = ...
    1 + s4;


d4 = ...
    s4*X_star/r4;



%% =========================================================
% Community matrix
%% =========================================================

M4 = ...
    r4/(1+s4*X_star)^2;



%% =========================================================
% Dominant eigenvalue
%% =========================================================

eig4 = ...
    eig(M4);


rho4 = ...
    max(abs(eig4));



%% =========================================================
% Stability
%% =========================================================

stab4 = ...
    -log(rho4);



%% =========================================================
% Simulation
%% =========================================================

X4 = ...
    zeros(1,max_step+1);


X4(1) = ...
    X0;



for ii = 1:max_step


    if X4(ii) == 0


        X4(ii+1) = 0;


    else


        X4(ii+1) = ...
            r4*X4(ii) ...
            /(1+s4*X4(ii));


        if X4(ii+1) < extinction_threshold


            X4(ii+1) = 0;


        end


    end


end



%% =========================================================
% Distance
%% =========================================================

D4 = ...
    abs(X4-X_star);



%% =========================================================
% Recovery time
%% =========================================================

recovery4 = ...
    NaN;


for ii = 1:length(X4)


    if D4(ii) <= recovery_threshold


        recovery4 = ...
            k(ii);


        break


    end


end



%% =========================================================
% Print parameters and recovery times
%% =========================================================

fprintf('\n');

fprintf('=============================================\n');

fprintf('Single-species Beverton-Holt model\n');

fprintf(...
    'X* = %.2f, X0 = %.2f\n',...
    X_star,...
    X0);


fprintf(...
    'Initial perturbation D0 = %.6f\n',...
    D0);


fprintf(...
    'Recovery threshold = %.6f\n',...
    recovery_threshold);


fprintf(...
    'Recovery criterion = 0.5%% of initial perturbation\n');


fprintf('=============================================\n');


fprintf(...
    ['s = %.2f, ',...
     'r = %.4f, ',...
     'd = %.4f, ',...
     'rho = %.4f, ',...
     'stability = %.4f, ',...
     'recovery time = %.0f\n'],...
     s1,...
     r1,...
     d1,...
     rho1,...
     stab1,...
     recovery1);


fprintf(...
    ['s = %.2f, ',...
     'r = %.4f, ',...
     'd = %.4f, ',...
     'rho = %.4f, ',...
     'stability = %.4f, ',...
     'recovery time = %.0f\n'],...
     s2,...
     r2,...
     d2,...
     rho2,...
     stab2,...
     recovery2);


fprintf(...
    ['s = %.2f, ',...
     'r = %.4f, ',...
     'd = %.4f, ',...
     'rho = %.4f, ',...
     'stability = %.4f, ',...
     'recovery time = %.0f\n'],...
     s3,...
     r3,...
     d3,...
     rho3,...
     stab3,...
     recovery3);


fprintf(...
    ['s = %.2f, ',...
     'r = %.4f, ',...
     'd = %.4f, ',...
     'rho = %.4f, ',...
     'stability = %.4f, ',...
     'recovery time = %.0f\n'],...
     s4,...
     r4,...
     d4,...
     rho4,...
     stab4,...
     recovery4);



%% =========================================================
%% Figure
%
% First row:
%
%       4 abundance panels
%
% Each:
%
%       4.5 cm x 3 cm
%
%
% Second row:
%
%       recovery time + stability
%
%       4.5 cm x 3 cm
%% =========================================================

figure


set(...
    gcf,...
    'unit','centimeters',...
    'position',[3 2 28 10.5]);



%% =========================================================
%% subplot 1
%% =========================================================

subplot(2,4,1)


set(...
    gca,...
    'unit','centimeters',...
    'position',[2 6.2 4.5 3]);


hold on



%% ---------------------------------------------------------
% Population trajectory
%% ---------------------------------------------------------

plot(...
    k,...
    X1,...
    '-o',...
    'Color',[46 167 224]/255,...
    'LineWidth',1,...
    'MarkerSize',3);



%% ---------------------------------------------------------
% Equilibrium
%% ---------------------------------------------------------

yline(...
    X_star,...
    ':',...
    'LineWidth',0.8);



%% ---------------------------------------------------------
% Recovery time
%% ---------------------------------------------------------

if ~isnan(recovery1)


    plot(...
        [recovery1 recovery1],...
        [0 X1(recovery1+1)],...
        '--',...
        'Color',[0.35 0.35 0.35],...
        'LineWidth',0.7);


    text(...
        recovery1/2,...
        0.10,...
        ['Recovery time = ' num2str(recovery1)],...
        'HorizontalAlignment','center',...
        'FontSize',8);


end



%% ---------------------------------------------------------
% Stability
%% ---------------------------------------------------------

text(...
    0.96,...
    0.88,...
    ['Stability = ' num2str(stab1,'%.3f')],...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'FontSize',8);



xlabel(...
    'Time step, $k$',...
    'Interpreter','latex');


ylabel(...
    'Abundance');



title(...
    ['$s^{\rm BH}=' num2str(s1,'%.2f')...
    ',\ d=' num2str(d1,'%.3f') '$'],...
    'Interpreter','latex',...
    'FontWeight','normal');



xlim([0 20]);

ylim([0 2.05]);


set(...
    gca,...
    'FontSize',8,...
    'Box','on');



%% =========================================================
%% subplot 2
%% =========================================================

subplot(2,4,2)


set(...
    gca,...
    'unit','centimeters',...
    'position',[8.5 6.2 4.5 3]);


hold on



plot(...
    k,...
    X2,...
    '-o',...
    'Color',[46 167 224]/255,...
    'LineWidth',1,...
    'MarkerSize',3);



yline(...
    X_star,...
    ':',...
    'LineWidth',0.8);



if ~isnan(recovery2)


    plot(...
        [recovery2 recovery2],...
        [0 X2(recovery2+1)],...
        '--',...
        'Color',[0.35 0.35 0.35],...
        'LineWidth',0.7);


    text(...
        recovery2/2,...
        0.10,...
        ['Recovery time = ' num2str(recovery2)],...
        'HorizontalAlignment','center',...
        'FontSize',8);


end



text(...
    0.96,...
    0.88,...
    ['Stability = ' num2str(stab2,'%.3f')],...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'FontSize',8);



xlabel(...
    'Time step, $k$',...
    'Interpreter','latex');


ylabel(...
    'Abundance');



title(...
    ['$s^{\rm BH}=' num2str(s2,'%.2f')...
    ',\ d=' num2str(d2,'%.3f') '$'],...
    'Interpreter','latex',...
    'FontWeight','normal');



xlim([0 20]);

ylim([0 2.05]);


set(...
    gca,...
    'FontSize',8,...
    'Box','on');



%% =========================================================
%% subplot 3
%% =========================================================

subplot(2,4,3)


set(...
    gca,...
    'unit','centimeters',...
    'position',[15 6.2 4.5 3]);


hold on



plot(...
    k,...
    X3,...
    '-o',...
    'Color',[46 167 224]/255,...
    'LineWidth',1,...
    'MarkerSize',3);



yline(...
    X_star,...
    ':',...
    'LineWidth',0.8);



if ~isnan(recovery3)


    plot(...
        [recovery3 recovery3],...
        [0 X3(recovery3+1)],...
        '--',...
        'Color',[0.35 0.35 0.35],...
        'LineWidth',0.7);


    text(...
        recovery3/2,...
        0.10,...
        ['Recovery time = ' num2str(recovery3)],...
        'HorizontalAlignment','center',...
        'FontSize',8);


end



text(...
    0.96,...
    0.88,...
    ['Stability = ' num2str(stab3,'%.3f')],...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'FontSize',8);



xlabel(...
    'Time step, $k$',...
    'Interpreter','latex');


ylabel(...
    'Abundance');



title(...
    ['$s^{\rm BH}=' num2str(s3,'%.2f')...
    ',\ d=' num2str(d3,'%.3f') '$'],...
    'Interpreter','latex',...
    'FontWeight','normal');



xlim([0 20]);

ylim([0 2.05]);


set(...
    gca,...
    'FontSize',8,...
    'Box','on');



%% =========================================================
%% subplot 4
%% =========================================================

subplot(2,4,4)


set(...
    gca,...
    'unit','centimeters',...
    'position',[21.5 6.2 4.5 3]);


hold on



plot(...
    k,...
    X4,...
    '-o',...
    'Color',[46 167 224]/255,...
    'LineWidth',1,...
    'MarkerSize',3);



yline(...
    X_star,...
    ':',...
    'LineWidth',0.8);



if ~isnan(recovery4)


    plot(...
        [recovery4 recovery4],...
        [0 X4(recovery4+1)],...
        '--',...
        'Color',[0.35 0.35 0.35],...
        'LineWidth',0.7);


    text(...
        recovery4/2,...
        0.10,...
        ['Recovery time = ' num2str(recovery4)],...
        'HorizontalAlignment','center',...
        'FontSize',8);


end



text(...
    0.96,...
    0.88,...
    ['Stability = ' num2str(stab4,'%.3f')],...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'FontSize',8);



xlabel(...
    'Time step, $k$',...
    'Interpreter','latex');


ylabel(...
    'Abundance');



title(...
    ['$s^{\rm BH}=' num2str(s4,'%.2f')...
    ',\ d=' num2str(d4,'%.3f') '$'],...
    'Interpreter','latex',...
    'FontWeight','normal');



xlim([0 20]);

ylim([0 2.05]);


set(...
    gca,...
    'FontSize',8,...
    'Box','on');



%% =========================================================
%% subplot 5
%
% Recovery time + stability
%
% Same format as the multi-species BH model.
%% =========================================================

subplot(2,4,5)


set(...
    gca,...
    'unit','centimeters',...
    'position',[11.75 1.2 4.5 3]);


hold on



%% =========================================================
% Prepare data
%% =========================================================

x_position = ...
    1:4;


recovery_time = ...
    [recovery1 recovery2 recovery3 recovery4];


stability_alpha = ...
    [stab1 stab2 stab3 stab4];


d_values = ...
    [d1 d2 d3 d4];



%% =========================================================
% LEFT AXIS
%
% Recovery time
%% =========================================================

yyaxis left


b1 = ...
    bar(...
    x_position-0.14,...
    recovery_time,...
    0.25);



ylabel(...
    'Recovery time',...
    'FontSize',9);



%% =========================================================
% Left y-axis range
%% =========================================================

finite_recovery = ...
    recovery_time(...
    ~isnan(recovery_time));



if isempty(finite_recovery)


    recovery_ymax = ...
        10;


else


    recovery_ymax = ...
        max(finite_recovery);


    recovery_ymax = ...
        max(recovery_ymax,1);


end



ylim(...
    [0 1.25*recovery_ymax]);



%% =========================================================
% Mark NR if necessary
%% =========================================================

for ii = 1:4


    if isnan(recovery_time(ii))


        text(...
            x_position(ii)-0.14,...
            1.08*recovery_ymax,...
            'NR',...
            'HorizontalAlignment','center',...
            'FontSize',7);


    end


end



%% =========================================================
% RIGHT AXIS
%
% Stability
%% =========================================================

yyaxis right


b2 = ...
    bar(...
    x_position+0.14,...
    stability_alpha,...
    0.25);



ylabel(...
    'Stability',...
    'FontSize',9);



%% =========================================================
% Stability boundary
%% =========================================================

yline(...
    0,...
    ':',...
    'LineWidth',0.5);



%% =========================================================
% Right y-axis range
%% =========================================================

stability_ymax = ...
    max(stability_alpha);


if stability_ymax <= 0


    stability_ymax = 1;


end


ylim(...
    [0 1.15*stability_ymax]);



%% =========================================================
% x axis
%% =========================================================

d_labels = ...
    {num2str(d1,'%.3f'),...
     num2str(d2,'%.3f'),...
     num2str(d3,'%.3f'),...
     num2str(d4,'%.3f')};



set(...
    gca,...
    'XTick',...
    x_position,...
    'XTickLabel',...
    d_labels,...
    'FontSize',8);



xlabel(...
    '$d$',...
    'Interpreter','latex',...
    'FontSize',9);



xlim([0.5 4.5])



%% =========================================================
% Legend
%% =========================================================

legend(...
    [b1 b2],...
    {'Recovery time','Stability'},...
    'Location','best',...
    'Box','off',...
    'FontSize',6);

title('Single-species BH model','FontSize',9); 

%% =========================================================
% Force x axis to black
%% =========================================================

ax = ...
    gca;


ax.XAxis.Color = ...
    [0 0 0];


ax.XAxis.Label.Color = ...
    [0 0 0];


box on