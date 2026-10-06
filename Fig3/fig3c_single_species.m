
% This file can be used to generate Fig. 3c

%% =========================================================
% Single-species Ricker model
%
% Recovery criterion:
%
% D(t) = |x(t) - x*|
%
% Recovery occurs when
%
% D(t) <= 0.005 * D0
%
% and the population remains within this range afterwards.
%
% Here:
% x* = 1
%
% =========================================================

clear all
clc
close all


%% =========================================================
% Common parameters
%% =========================================================

x0 = 1.0;             % equilibrium abundance x*

T = 20;
dt = 1;

steps = T/dt;

time = 0:dt:T;


% ---------------------------------------------------------
% disturbance_index = 4 means:
%
% x(4) corresponds to time = 3
%
% so the disturbance is introduced at t = 3
% ---------------------------------------------------------

disturbance_index = 4;

disturbance_value = 0.5;


% ---------------------------------------------------------
% Recovery criterion
%
% 0.5% = 0.005
% ---------------------------------------------------------

recovery_ratio = 0.005;


%% =========================================================
% s range and theoretical stability
%% =========================================================

s_all_1 = 0:0.01:0.99;

s_all_2 = 1.01:0.01:4.2;

s_all = [s_all_1, s_all_2];

stability_all = -log(abs(1-s_all));


%% =========================================================
% Figure
%% =========================================================

figure

set(gcf,...
    'unit','centimeters',...
    'position',[10 5 18 15]);



%% =========================================================
%% subplot 1
%
% Low self-regulation
%
% s = 0.4
%% =========================================================

subplot(2,2,1)

set(gca,...
    'unit','centimeters',...
    'position',[2 8 4.5 3]);


s = 0.4;

r = s;


%% ---------------------------------------------------------
% Simulation
%% ---------------------------------------------------------

x1 = zeros(1,steps+1);

x1(1) = x0;


for k = 1:steps

    if k == disturbance_index

        x1(k) = disturbance_value;

    end


    x1(k+1) = ...
        x1(k)*exp(r-s*x1(k));

end


%% ---------------------------------------------------------
% Stability
%
% lambda = 1-s
%
% alpha = -log(|lambda|)
%% ---------------------------------------------------------

stab1 = -log(abs(1-s));


%% ---------------------------------------------------------
% Recovery time
%% ---------------------------------------------------------

distance1 = abs(x1-x0);


% initial disturbance magnitude
initial_disturbance1 = ...
    distance1(disturbance_index);


% recovery threshold
threshold1 = ...
    recovery_ratio*initial_disturbance1;


recovery_step1 = NaN;

recovery_time1 = NaN;


for kk = disturbance_index+1:length(x1)

    % First time entering the recovery region
    % AND remaining in the region afterwards

    if all(distance1(kk:end) <= threshold1)

        recovery_step1 = kk;

        recovery_time1 = ...
            time(kk)-time(disturbance_index);

        break

    end

end


fprintf('\n==============================\n')
fprintf('s = %.2f\n',s)
fprintf('Stability = %.4f\n',stab1)
fprintf('Initial disturbance = %.4f\n',initial_disturbance1)
fprintf('Recovery threshold = %.6f\n',threshold1)

if isnan(recovery_time1)

    fprintf('Population did not recover.\n')

else

    fprintf('Recovery time = %.0f generations\n',...
        recovery_time1)

    fprintf('Recovery occurs at t = %.0f\n',...
        time(recovery_step1))

end



%% ---------------------------------------------------------
% Plot
%% ---------------------------------------------------------

hold on


plot(...
    time,...
    x1,...
    '-o',...
    'Color',[46 167 224]/255,...
    'LineWidth',1,...
    'MarkerSize',3);


grid off

box on


xlim([0 20])

ylim([-0.1 2.1])


yticks(0:0.5:2.5)

yticklabels(...
    {'0' '' '1' '' '2' ''});


xticks(0:1:20)

xticklabels(...
    {'0' '' '' '' '' '5' '' '' '' '' ...
     '10' '' '' '' '' '15' '' '' '' '' '20'});



% disturbance
line(...
    [3 3],...
    [-0.1 0.5],...
    'LineStyle','--',...
    'LineWidth',0.5,...
    'Color','k');


ylabel(...
    'Abundance',...
    'fontsize',9);


title(...
    'Low level of self-regulation',...
    'fontsize',9);



%% =========================================================
%% subplot 2
%
% Moderate self-regulation
%
% s = 0.9
%% =========================================================

subplot(2,2,2)

set(gca,...
    'unit','centimeters',...
    'position',[9 8 4.5 3]);


s = 0.9;

r = s;


%% ---------------------------------------------------------
% Simulation
%% ---------------------------------------------------------

x2 = zeros(1,steps+1);

x2(1) = x0;


for k = 1:steps

    if k == disturbance_index

        x2(k) = disturbance_value;

    end


    x2(k+1) = ...
        x2(k)*exp(r-s*x2(k));

end


%% ---------------------------------------------------------
% Stability
%% ---------------------------------------------------------

stab2 = -log(abs(1-s));


%% ---------------------------------------------------------
% Recovery time
%% ---------------------------------------------------------

distance2 = abs(x2-x0);


initial_disturbance2 = ...
    distance2(disturbance_index);


threshold2 = ...
    recovery_ratio*initial_disturbance2;


recovery_step2 = NaN;

recovery_time2 = NaN;


for kk = disturbance_index+1:length(x2)

    if all(distance2(kk:end) <= threshold2)

        recovery_step2 = kk;

        recovery_time2 = ...
            time(kk)-time(disturbance_index);

        break

    end

end


fprintf('\n==============================\n')
fprintf('s = %.2f\n',s)
fprintf('Stability = %.4f\n',stab2)
fprintf('Initial disturbance = %.4f\n',initial_disturbance2)
fprintf('Recovery threshold = %.6f\n',threshold2)

if isnan(recovery_time2)

    fprintf('Population did not recover.\n')

else

    fprintf('Recovery time = %.0f generations\n',...
        recovery_time2)

    fprintf('Recovery occurs at t = %.0f\n',...
        time(recovery_step2))

end



%% ---------------------------------------------------------
% Plot
%% ---------------------------------------------------------

hold on


plot(...
    time,...
    x2,...
    '-o',...
    'Color',[46 167 224]/255,...
    'LineWidth',1,...
    'MarkerSize',3);


grid off

box on


xlim([0 20])

ylim([-0.1 2.1])


yticks(0:0.5:2.5)

yticklabels(...
    {'0' '' '1' '' '2' ''});


xticks(0:1:20)

xticklabels(...
    {'0' '' '' '' '' '5' '' '' '' '' ...
     '10' '' '' '' '' '15' '' '' '' '' '20'});



line(...
    [3 3],...
    [-0.1 0.5],...
    'LineStyle','--',...
    'LineWidth',0.5,...
    'Color','k');


title(...
    'Moderate level of self-regulation',...
    'fontsize',9);



%% =========================================================
%% subplot 3
%
% High self-regulation
%
% s = 1.73
%% =========================================================

subplot(2,2,3)

set(gca,...
    'unit','centimeters',...
    'position',[2 3 4.5 3]);


s = 1.73;

r = s;


%% ---------------------------------------------------------
% Simulation
%% ---------------------------------------------------------

x3 = zeros(1,steps+1);

x3(1) = x0;


for k = 1:steps

    if k == disturbance_index

        x3(k) = disturbance_value;

    end


    x3(k+1) = ...
        x3(k)*exp(r-s*x3(k));

end


%% ---------------------------------------------------------
% Stability
%% ---------------------------------------------------------

stab3 = -log(abs(1-s));


%% ---------------------------------------------------------
% Recovery time
%% ---------------------------------------------------------

distance3 = abs(x3-x0);


initial_disturbance3 = ...
    distance3(disturbance_index);


threshold3 = ...
    recovery_ratio*initial_disturbance3;


recovery_step3 = NaN;

recovery_time3 = NaN;


for kk = disturbance_index+1:length(x3)

    if all(distance3(kk:end) <= threshold3)

        recovery_step3 = kk;

        recovery_time3 = ...
            time(kk)-time(disturbance_index);

        break

    end

end


fprintf('\n==============================\n')
fprintf('s = %.2f\n',s)
fprintf('Stability = %.4f\n',stab3)
fprintf('Initial disturbance = %.4f\n',initial_disturbance3)
fprintf('Recovery threshold = %.6f\n',threshold3)

if isnan(recovery_time3)

    fprintf('Population did not recover.\n')

else

    fprintf('Recovery time = %.0f generations\n',...
        recovery_time3)

    fprintf('Recovery occurs at t = %.0f\n',...
        time(recovery_step3))

end



%% ---------------------------------------------------------
% Plot
%% ---------------------------------------------------------

hold on


plot(...
    time,...
    x3,...
    '-o',...
    'Color',[46 167 224]/255,...
    'LineWidth',1,...
    'MarkerSize',3);


grid off

box on


xlim([0 20])

ylim([-0.1 2.1])


yticks(0:0.5:2.5)

yticklabels(...
    {'0' '' '1' '' '2' ''});


xticks(0:1:20)

xticklabels(...
    {'0' '' '' '' '' '5' '' '' '' '' ...
     '10' '' '' '' '' '15' '' '' '' '' '20'});



line(...
    [3 3],...
    [-0.1 0.5],...
    'LineStyle','--',...
    'LineWidth',0.5,...
    'Color','k');


ylabel(...
    'Abundance',...
    'fontsize',9);


title(...
    'High level of self-regulation',...
    'fontsize',9);



%% =========================================================
%% subplot 4
%
% Extremely high self-regulation
%
% s = 4
%% =========================================================

subplot(2,2,4)

set(gca,...
    'unit','centimeters',...
    'position',[9 3 4.5 3]);


s = 4.0;

r = s;


%% ---------------------------------------------------------
% Simulation
%% ---------------------------------------------------------

x4 = zeros(1,steps+1);

x4(1) = x0;


for k = 1:steps

    if k == disturbance_index

        x4(k) = disturbance_value;

    end


    x4(k+1) = ...
        x4(k)*exp(r-s*x4(k));


    % extinction criterion
    if x4(k+1) < 1e-4

        x4(k+1:end) = 0;

        break

    end

end


%% ---------------------------------------------------------
% Stability
%% ---------------------------------------------------------

stab4 = -log(abs(1-s));


%% ---------------------------------------------------------
% Recovery time
%% ---------------------------------------------------------

distance4 = abs(x4-x0);


initial_disturbance4 = ...
    distance4(disturbance_index);


threshold4 = ...
    recovery_ratio*initial_disturbance4;


recovery_step4 = NaN;

recovery_time4 = NaN;


for kk = disturbance_index+1:length(x4)

    if all(distance4(kk:end) <= threshold4)

        recovery_step4 = kk;

        recovery_time4 = ...
            time(kk)-time(disturbance_index);

        break

    end

end


fprintf('\n==============================\n')
fprintf('s = %.2f\n',s)
fprintf('Stability = %.4f\n',stab4)
fprintf('Initial disturbance = %.4f\n',initial_disturbance4)
fprintf('Recovery threshold = %.6f\n',threshold4)

if isnan(recovery_time4)

    fprintf('Population did not recover.\n')

else

    fprintf('Recovery time = %.0f generations\n',...
        recovery_time4)

    fprintf('Recovery occurs at t = %.0f\n',...
        time(recovery_step4))

end



%% ---------------------------------------------------------
% Plot
%% ---------------------------------------------------------

hold on


plot(...
    time,...
    x4,...
    '-o',...
    'Color',[46 167 224]/255,...
    'LineWidth',1,...
    'MarkerSize',3);


grid off

box on

xlim([0 20])

ylim([-0.1 4.1])


yticks(0:0.5:4)

yticklabels(...
    {'0' '' '1' '' '2' '' '3' '' '4'});


xticks(0:1:20)

xticklabels(...
    {'0' '' '' '' '' '5' '' '' '' '' ...
     '10' '' '' '' '' '15' '' '' '' '' '20'});


line(...
    [3 3],...
    [-0.1 0.5],...
    'LineStyle','--',...
    'LineWidth',0.5,...
    'Color','k');


title(...
    'Extremely high level of self-regulation',...
    'fontsize',9);