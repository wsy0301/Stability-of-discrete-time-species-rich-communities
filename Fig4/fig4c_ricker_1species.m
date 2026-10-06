
% This file can be used to generate Fig. 4c 


%% =========================================================
% Single-species Ricker model
%
% Model:
%
%       X(k+1)
%       =
%       X(k)*exp[r-sX(k)]
%
%
% Prescribed equilibrium:
%
%       X* = 1
%
%
% Equilibrium condition:
%
%       r-sX* = 0
%
% Since X*=1:
%
%       r = s
%
%
% Self-regulation:
%
%       d = sX*
%
% Since X*=1:
%
%       d = s
%
%
% Therefore:
%
%       r = s = d
%
%
% Jacobian / community matrix:
%
%       M
%       =
%       d/dX [X exp(r-sX)] at X=X*
%
%       =
%       exp(r-sX*)[1-sX*]
%
% At equilibrium:
%
%       exp(r-sX*) = 1
%
% Therefore:
%
%       M = 1-sX*
%
% Since X*=1 and d=s:
%
%       M = 1-d
%
%
% Spectral radius:
%
%       rho = |M|
%           = |1-d|
%
%
% Stability:
%
%       alpha = -log(rho)
%
%
% Recovery criterion:
%
%       D0 = |X(0)-X*|
%
%       D(k) = |X(k)-X*|
%
% Recovered when:
%
%       D(k) <= 0.005 D0
%
% i.e. perturbation has declined to
% 0.5% of its initial magnitude.
%
%% =========================================================

clear all
clc
close all



%% =========================================================
% Common parameters
%% =========================================================

% Equilibrium abundance
X_star = 1;


% Initial abundance
X0 = 0.8;


% Maximum simulation time
max_step = 100;


% Time vector
k = ...
    0:max_step;


% Extinction threshold
extinction_threshold = ...
    0.0001;



%% =========================================================
% Recovery criterion
%% =========================================================

% Initial perturbation magnitude
D0 = ...
    abs(X0-X_star);


% 0.5% of initial perturbation
recovery_fraction = ...
    0.005;


% Recovery threshold
recovery_threshold = ...
    recovery_fraction*D0;



fprintf('\n')

disp('============================================')
disp('RECOVERY CRITERION')
disp('============================================')


fprintf(...
    'X* = %.4f\n',...
    X_star);


fprintf(...
    'X0 = %.4f\n',...
    X0);


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
% Four self-regulation strengths
%
% Same values as the original Ricker program.
%% =========================================================

d1 = 0.5;

d2 = 0.8;

d3 = 1.4;

d4 = 1.8;



%% =========================================================
%% CASE 1
%
% d = 0.5
%% =========================================================

% Since X*=1:
%
%       s = d
%       r = d

s1 = ...
    d1;


r1 = ...
    d1;



%% ---------------------------------------------------------
% Community matrix
%
%       M = 1-d
%% ---------------------------------------------------------

M1 = ...
    1-d1;



%% ---------------------------------------------------------
% Spectral radius
%% ---------------------------------------------------------

rho1 = ...
    abs(M1);



%% ---------------------------------------------------------
% Stability
%% ---------------------------------------------------------

if rho1 == 0


    alpha1 = ...
        Inf;


else


    alpha1 = ...
        -log(rho1);


end



%% ---------------------------------------------------------
% Simulation
%% ---------------------------------------------------------

X1 = ...
    zeros(1,max_step+1);


X1(1) = ...
    X0;



for ii = 1:max_step


    if X1(ii) == 0


        X1(ii+1) = 0;


    else


        X1(ii+1) = ...
            X1(ii) ...
            *exp(...
            r1 ...
            - s1*X1(ii));


        if X1(ii+1) < extinction_threshold


            X1(ii+1) = 0;


        end


    end


end



%% ---------------------------------------------------------
% Distance from equilibrium
%% ---------------------------------------------------------

D1 = ...
    abs(X1-X_star);



%% ---------------------------------------------------------
% Recovery time
%% ---------------------------------------------------------

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
%% CASE 2
%
% d = 0.8
%% =========================================================

s2 = ...
    d2;


r2 = ...
    d2;



%% ---------------------------------------------------------
% Community matrix
%% ---------------------------------------------------------

M2 = ...
    1-d2;



%% ---------------------------------------------------------
% Spectral radius
%% ---------------------------------------------------------

rho2 = ...
    abs(M2);



%% ---------------------------------------------------------
% Stability
%% ---------------------------------------------------------

if rho2 == 0


    alpha2 = ...
        Inf;


else


    alpha2 = ...
        -log(rho2);


end



%% ---------------------------------------------------------
% Simulation
%% ---------------------------------------------------------

X2 = ...
    zeros(1,max_step+1);


X2(1) = ...
    X0;



for ii = 1:max_step


    if X2(ii) == 0


        X2(ii+1) = 0;


    else


        X2(ii+1) = ...
            X2(ii) ...
            *exp(...
            r2 ...
            - s2*X2(ii));


        if X2(ii+1) < extinction_threshold


            X2(ii+1) = 0;


        end


    end


end



%% ---------------------------------------------------------
% Distance
%% ---------------------------------------------------------

D2 = ...
    abs(X2-X_star);



%% ---------------------------------------------------------
% Recovery time
%% ---------------------------------------------------------

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
%% CASE 3
%
% d = 1.4
%% =========================================================

s3 = ...
    d3;


r3 = ...
    d3;



%% ---------------------------------------------------------
% Community matrix
%% ---------------------------------------------------------

M3 = ...
    1-d3;



%% ---------------------------------------------------------
% Spectral radius
%% ---------------------------------------------------------

rho3 = ...
    abs(M3);



%% ---------------------------------------------------------
% Stability
%% ---------------------------------------------------------

if rho3 == 0


    alpha3 = ...
        Inf;


else


    alpha3 = ...
        -log(rho3);


end



%% ---------------------------------------------------------
% Simulation
%% ---------------------------------------------------------

X3 = ...
    zeros(1,max_step+1);


X3(1) = ...
    X0;



for ii = 1:max_step


    if X3(ii) == 0


        X3(ii+1) = 0;


    else


        X3(ii+1) = ...
            X3(ii) ...
            *exp(...
            r3 ...
            - s3*X3(ii));


        if X3(ii+1) < extinction_threshold


            X3(ii+1) = 0;


        end


    end


end



%% ---------------------------------------------------------
% Distance
%% ---------------------------------------------------------

D3 = ...
    abs(X3-X_star);



%% ---------------------------------------------------------
% Recovery time
%% ---------------------------------------------------------

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
%% CASE 4
%
% d = 1.8
%% =========================================================

s4 = ...
    d4;


r4 = ...
    d4;



%% ---------------------------------------------------------
% Community matrix
%% ---------------------------------------------------------

M4 = ...
    1-d4;



%% ---------------------------------------------------------
% Spectral radius
%% ---------------------------------------------------------

rho4 = ...
    abs(M4);



%% ---------------------------------------------------------
% Stability
%% ---------------------------------------------------------

if rho4 == 0


    alpha4 = ...
        Inf;


else


    alpha4 = ...
        -log(rho4);


end



%% ---------------------------------------------------------
% Simulation
%% ---------------------------------------------------------

X4 = ...
    zeros(1,max_step+1);


X4(1) = ...
    X0;



for ii = 1:max_step


    if X4(ii) == 0


        X4(ii+1) = 0;


    else


        X4(ii+1) = ...
            X4(ii) ...
            *exp(...
            r4 ...
            - s4*X4(ii));


        if X4(ii+1) < extinction_threshold


            X4(ii+1) = 0;


        end


    end


end



%% ---------------------------------------------------------
% Distance
%% ---------------------------------------------------------

D4 = ...
    abs(X4-X_star);



%% ---------------------------------------------------------
% Recovery time
%% ---------------------------------------------------------

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
% Print results
%% =========================================================

fprintf('\n')

disp('============================================')
disp('SINGLE-SPECIES RICKER MODEL')
disp('============================================')


fprintf(...
    ['d = %.2f, ',...
     'r = %.2f, ',...
     's = %.2f, ',...
     'M = %.4f, ',...
     'rho = %.4f, ',...
     'alpha = %.4f, ',...
     'recovery time = %.0f\n'],...
     d1,...
     r1,...
     s1,...
     M1,...
     rho1,...
     alpha1,...
     recovery1);



fprintf(...
    ['d = %.2f, ',...
     'r = %.2f, ',...
     's = %.2f, ',...
     'M = %.4f, ',...
     'rho = %.4f, ',...
     'alpha = %.4f, ',...
     'recovery time = %.0f\n'],...
     d2,...
     r2,...
     s2,...
     M2,...
     rho2,...
     alpha2,...
     recovery2);



fprintf(...
    ['d = %.2f, ',...
     'r = %.2f, ',...
     's = %.2f, ',...
     'M = %.4f, ',...
     'rho = %.4f, ',...
     'alpha = %.4f, ',...
     'recovery time = %.0f\n'],...
     d3,...
     r3,...
     s3,...
     M3,...
     rho3,...
     alpha3,...
     recovery3);



fprintf(...
    ['d = %.2f, ',...
     'r = %.2f, ',...
     's = %.2f, ',...
     'M = %.4f, ',...
     'rho = %.4f, ',...
     'alpha = %.4f, ',...
     'recovery time = %.0f\n'],...
     d4,...
     r4,...
     s4,...
     M4,...
     rho4,...
     alpha4,...
     recovery4);



%% =========================================================
% Values used for plotting
%% =========================================================

Tshow = ...
    30;


time_show = ...
    0:Tshow;



%% =========================================================
%% ONE FIGURE
%
% First row:
%
%       four abundance panels
%
% Second row:
%
%       recovery time + stability
%
% Every panel:
%
%       width  = 4.5 cm
%       height = 3 cm
%% =========================================================

figure


set(...
    gcf,...
    'unit','centimeters',...
    'position',[3 2 28 10.5]);



%% =========================================================
%% subplot 1
%
% d = 0.5
%% =========================================================

subplot(2,4,1)


set(...
    gca,...
    'unit','centimeters',...
    'position',[2 6.2 4.5 3]);


hold on



plot(...
    time_show,...
    X1(1:Tshow+1),...
    '-o',...
    'Color',[46 167 224]/255,...
    'LineWidth',1,...
    'MarkerSize',3);



% Equilibrium
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


end



text(...
    0.96,...
    0.08,...
    sprintf(...
    '$d=%.1f,\\;T_R=%.0f,\\;\\alpha=%.3f$',...
    d1,...
    recovery1,...
    alpha1),...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Interpreter','latex',...
    'FontSize',7);



xlabel(...
    'Time step, $k$',...
    'Interpreter','latex');


ylabel(...
    'Abundance');



xlim([0 Tshow])

ylim([0.75 1.18])


set(...
    gca,...
    'FontSize',8,...
    'Box','on');



%% =========================================================
%% subplot 2
%
% d = 0.8
%% =========================================================

subplot(2,4,2)


set(...
    gca,...
    'unit','centimeters',...
    'position',[8.5 6.2 4.5 3]);


hold on



plot(...
    time_show,...
    X2(1:Tshow+1),...
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


end



text(...
    0.96,...
    0.08,...
    sprintf(...
    '$d=%.1f,\\;T_R=%.0f,\\;\\alpha=%.3f$',...
    d2,...
    recovery2,...
    alpha2),...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Interpreter','latex',...
    'FontSize',7);



xlabel(...
    'Time step, $k$',...
    'Interpreter','latex');


ylabel(...
    'Abundance');



xlim([0 Tshow])

ylim([0.75 1.18])


set(...
    gca,...
    'FontSize',8,...
    'Box','on');



%% =========================================================
%% subplot 3
%
% d = 1.4
%% =========================================================

subplot(2,4,3)


set(...
    gca,...
    'unit','centimeters',...
    'position',[15 6.2 4.5 3]);


hold on



plot(...
    time_show,...
    X3(1:Tshow+1),...
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


end



text(...
    0.96,...
    0.08,...
    sprintf(...
    '$d=%.1f,\\;T_R=%.0f,\\;\\alpha=%.3f$',...
    d3,...
    recovery3,...
    alpha3),...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Interpreter','latex',...
    'FontSize',7);



xlabel(...
    'Time step, $k$',...
    'Interpreter','latex');


ylabel(...
    'Abundance');



xlim([0 Tshow])

ylim([0.75 1.18])


set(...
    gca,...
    'FontSize',8,...
    'Box','on');



%% =========================================================
%% subplot 4
%
% d = 1.8
%% =========================================================

subplot(2,4,4)


set(...
    gca,...
    'unit','centimeters',...
    'position',[21.5 6.2 4.5 3]);


hold on



plot(...
    time_show,...
    X4(1:Tshow+1),...
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


end



text(...
    0.96,...
    0.08,...
    sprintf(...
    '$d=%.1f,\\;T_R=%.0f,\\;\\alpha=%.3f$',...
    d4,...
    recovery4,...
    alpha4),...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Interpreter','latex',...
    'FontSize',7);



xlabel(...
    'Time step, $k$',...
    'Interpreter','latex');


ylabel(...
    'Abundance');



xlim([0 Tshow])

ylim([0.75 1.18])


set(...
    gca,...
    'FontSize',8,...
    'Box','on');



%% =========================================================
%% subplot 5
%
% Recovery time + stability
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
    [alpha1 alpha2 alpha3 alpha4];


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
    [0 1.20*recovery_ymax]);



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



stability_ymax = ...
    max(stability_alpha);


if isfinite(stability_ymax)


    ylim(...
        [0 1.45*stability_ymax]);


end



%% =========================================================
% x axis
%% =========================================================

set(...
    gca,...
    'XTick',...
    x_position,...
    'XTickLabel',...
    {'0.5','0.8','1.4','1.8'},...
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



title(...
    'Single-species Ricker model',...
    'FontSize',9);



%% =========================================================
% Force x-axis to black
%% =========================================================

ax = ...
    gca;


ax.XAxis.Color = ...
    [0 0 0];


ax.XAxis.Label.Color = ...
    [0 0 0];


box on