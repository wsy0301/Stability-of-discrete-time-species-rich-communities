
% This file can be used to generate Figs. 4a-b



clear all
clc
close all


%% =========================================================
% Parameters for multi-species communities
%% =========================================================

S = 100;          % number of species
C = 0.10;         % connectance
sigma = 0.01;     % interaction strength SD

rng(1)            % fix random seed for reproducibility


%% =========================================================
% Generate Foff for a competitive community
%
% Competitive:
%
%       F_ij < 0
%       F_ji < 0
%
% Each direction is independently generated.
%% =========================================================

Foff = zeros(S,S);

for i = 2:S
    for j = 1:i-1

        p1 = rand;

        if p1 <= C

            Foff(i,j) = -abs(normrnd(0,sigma));
            Foff(j,i) = -abs(normrnd(0,sigma));

        end

    end
end

Foff(1:S+1:end) = 0;


%% =========================================================
% Row sums
%
% sumF_i = sum_{j~=i} F_ij
%% =========================================================

sumF = sum(Foff,2);


%% =========================================================
% 1. Single-species Ricker
%
% d = r
%% =========================================================

r_ricker_single = linspace(0,10,500);
d_ricker_single = r_ricker_single;


%% =========================================================
% 2. Multi-species Ricker
%
% r_i = d - sumF_i
%
% For each d, compute the proportion of species with r_i > 0
%% =========================================================

d_ricker_multi = linspace(0,10,500);
prop_ricker_multi = zeros(size(d_ricker_multi));

for k = 1:length(d_ricker_multi)

    d = d_ricker_multi(k);

    r_i = d - sumF;

    prop_ricker_multi(k) = sum(r_i > 0)/S;

end


%% =========================================================
% 3. Single-species Beverton-Holt
%
% d = (r-1)/r = 1 - 1/r
%
% Need r > 1
%% =========================================================

r_bh_single = linspace(1.01,20,500);
d_bh_single = (r_bh_single - 1)./r_bh_single;


%% =========================================================
% 4. Multi-species Beverton-Holt
%
% r_i = 1 / (1 - d + sumF_i)
%
% r_i > 0  <=>  1 - d + sumF_i > 0
%
% For each d, compute the proportion of species with r_i > 0
%% =========================================================

% Since sumF is negative for competitive community,
% positive-r region typically lies below d < 1 + sumF_i.
% We choose a d range that can show the transition clearly.

d_bh_multi = linspace(0,1.5,500);
prop_bh_multi = zeros(size(d_bh_multi));

for k = 1:length(d_bh_multi)

    d = d_bh_multi(k);

    denominator = 1 - d + sumF;

    prop_bh_multi(k) = sum(denominator > 0)/S;

end

%% =========================================================
% Maximum d that guarantees r_i > 0 for ALL species
%% =========================================================

dcrit_bh = min(1 + sumF);

fprintf('\n');
fprintf('For all species to have r_i > 0:\n');
fprintf('d must satisfy d < %.4f\n',dcrit_bh);

%% =========================================================
% Some useful diagnostic output
%% =========================================================

fprintf('\n');
disp('============================================');
disp('Competitive community diagnostics');
disp('============================================');

fprintf('S = %d\n',S);
fprintf('Target connectance C = %.3f\n',C);
fprintf('Realized connectance = %.3f\n',nnz(Foff)/(S*(S-1)));
fprintf('sigma = %.3f\n',sigma);

fprintf('\n');
fprintf('min(sumF) = %.4f\n',min(sumF));
fprintf('max(sumF) = %.4f\n',max(sumF));

fprintf('\n');
fprintf('For multi-species BH,\n');
fprintf('species i has positive r_i only if d < 1 + sumF_i.\n');
fprintf('min(1+sumF) = %.4f\n',min(1+sumF));
fprintf('max(1+sumF) = %.4f\n',max(1+sumF));


%% =========================================================
% Figure with 4 subplots
%
% Each subplot:
% width  = 4.5 cm
% height = 3 cm
%% =========================================================

figure
set(gcf,...
    'unit','centimeters',...
    'position',[3 3 28 8]);


%% =========================================================
% Subplot 1: single-species Ricker
%% =========================================================

subplot(1,4,1)
set(gca,...
    'unit','centimeters',...
    'position',[1 2 4.5 3]);
hold on

plot(r_ricker_single,...
     d_ricker_single,...
     'k-',...
     'LineWidth',1);

xlabel('$r$','Interpreter','latex','FontSize',9);
ylabel('$d$','Interpreter','latex','FontSize',9);

title('Single-species Ricker','FontSize',9);

xlim([0 5.1]);
ylim([0 5.1]);

set(gca,'FontSize',8);
box on


%% =========================================================
% Subplot 2: multi-species Ricker
%% =========================================================

subplot(1,4,2)
set(gca,...
    'unit','centimeters',...
    'position',[7 2 4.5 3]);
hold on

plot(d_ricker_multi,...
     prop_ricker_multi,...
     'k-',...
     'LineWidth',1);

xlabel('$d$','Interpreter','latex','FontSize',9);
ylabel('P($r_i>0$)','Interpreter','latex','FontSize',9);

title('Multi-species Ricker','FontSize',9);

xlim([0 5.1]);
ylim([0 1.05]);

set(gca,'FontSize',8);
box on


%% =========================================================
% Subplot 3: single-species Beverton-Holt
%% =========================================================

subplot(1,4,3)
set(gca,...
    'unit','centimeters',...
    'position',[13 2 4.5 3]);
hold on

plot(r_bh_single,...
     d_bh_single,...
     'k-',...
     'LineWidth',1);

xlabel('$r$','Interpreter','latex','FontSize',9);
ylabel('$d$','Interpreter','latex','FontSize',9);

title('Single-species BH','FontSize',9);

xlim([1 10.1]);
ylim([0 1]);

set(gca,'FontSize',8);
box on


%% =========================================================
% Subplot 4: multi-species Beverton-Holt
%% =========================================================

subplot(1,4,4)
set(gca,...
    'unit','centimeters',...
    'position',[19 2 4.5 3]);
hold on

plot(d_bh_multi,...
     prop_bh_multi,...
     'k-',...
     'LineWidth',1);
 
 xline(dcrit_bh,...
      '--',...
      ['$d_{\rm crit}=' num2str(dcrit_bh,'%.2f') '$'],...
      'Interpreter','latex',...
      'LabelVerticalAlignment','middle',...
      'LabelHorizontalAlignment','left',...
      'LineWidth',1);

xlabel('$d$','Interpreter','latex','FontSize',9);
ylabel('P($r_i>0$)','Interpreter','latex','FontSize',9);

title('Multi-species BH','FontSize',9);

xlim([0 1.4]);
ylim([0 1.05]);

set(gca,'FontSize',8);
box on