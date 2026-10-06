
% This file can be used to generate Fig. 7c



clear all
clc
close all


%% =========================================================
% Parameters
%% =========================================================

sigma = 0.1;
S = 200;
C = 0.1;
n = 100;

% Three self-regulation values
val_d = [0.5 1 1.5];

% Fix random seed
rng(1);


%% =========================================================
% Preallocate stability
%
% row 1: d = 0.5
% row 2: d = 1
% row 3: d = 1.5
%% =========================================================

J_cas   = zeros(3,n);
J_case1 = zeros(3,n);
J_case2 = zeros(3,n);
J_case3 = zeros(3,n);


%% =========================================================
%% 1. Cascade food web
%% =========================================================

for t = 1:n

    % -----------------------------------------------------
    % Generate off-diagonal community matrix
    % -----------------------------------------------------

    Moff = zeros(S);

    for i = 1:S
        for j = 1:S

            if i < j

                p1 = rand;

                if p1 <= C

                    Moff(i,j) = -abs(normrnd(0,sigma,1,1));
                    Moff(j,i) =  abs(normrnd(0,sigma,1,1));

                end

            end

        end
    end


    % -----------------------------------------------------
    % Eigenvalues of off-diagonal matrix
    % -----------------------------------------------------

    lambda_off = eig(Moff);


    % -----------------------------------------------------
    % Calculate stability for different d
    %
    % M = Moff + (1-d)I
    %
    % therefore
    %
    % lambda(M) = lambda(Moff) + (1-d)
    % -----------------------------------------------------

    for k = 1:length(val_d)

        d = val_d(k);

        lambda_M = lambda_off + (1-d);

        rho = max(abs(lambda_M));

        J_cas(k,t) = -log(rho);

    end

end


%% =========================================================
%% 2. Cascade + intervality
%% =========================================================

for t = 1:n


    % -----------------------------------------------------
    % Number of prey species
    % -----------------------------------------------------

    L = zeros(1,S);

    for i = 1:S

        L(i) = fix(C*(i-1));

    end


    % -----------------------------------------------------
    % Starting point of feeding interval
    % -----------------------------------------------------

    s = zeros(1,S);

    for i = 1:S

        s(i) = unidrnd(i-L(i));

    end


    % -----------------------------------------------------
    % Generate adjacency matrix K
    %
    % K(i,j) = 1:
    % species i preys on species j
    % -----------------------------------------------------

    K = zeros(S);

    for i = 1:S
        for j = 1:S

            if j >= s(i) && j <= s(i)+L(i)-1

                K(i,j) = 1;

            else

                K(i,j) = 0;

            end

        end
    end


    % -----------------------------------------------------
    % Community matrix
    % -----------------------------------------------------

    Moff = zeros(S);

    for i = 1:S
        for j = 1:S

            if K(i,j)==1 && K(j,i)==0

                Moff(i,j) =  abs(normrnd(0,sigma,1,1));
                Moff(j,i) = -abs(normrnd(0,sigma,1,1));

            end

        end
    end


    % -----------------------------------------------------
    % Eigenvalues
    % -----------------------------------------------------

    lambda_off = eig(Moff);


    % -----------------------------------------------------
    % Different d
    % -----------------------------------------------------

    for k = 1:length(val_d)

        d = val_d(k);

        lambda_M = lambda_off + (1-d);

        rho = max(abs(lambda_M));

        J_case1(k,t) = -log(rho);

    end

end


%% =========================================================
%% 3. Cascade + broad degree distribution
%% =========================================================

for t = 1:n


    % -----------------------------------------------------
    % Niche values eta
    % -----------------------------------------------------

    eta = zeros(1,S);

    for i = 1:S

        eta(i) = unifrnd(0,1);

    end

    eta = sort(eta);


    % -----------------------------------------------------
    % Niche radius
    % -----------------------------------------------------

    beta = zeros(1,S);
    r = zeros(1,S);

    for i = 1:S

        beta(i) = betarnd(1,1/C-1);

        r(i) = eta(i)*beta(i);

    end


    % -----------------------------------------------------
    % Number of prey species L
    % -----------------------------------------------------

    L = zeros(1,S);

    for i = 1:S

        L(i) = binornd(i-1,r(i));

    end


    % -----------------------------------------------------
    % Effective connectance
    % -----------------------------------------------------

    C_eff = zeros(1,S);

    C_eff(1) = 0;

    for i = 2:S

        C_eff(i) = L(i)/(i-1);

    end


    % -----------------------------------------------------
    % Generate adjacency matrix K
    % -----------------------------------------------------

    K = zeros(S);

    for i = 1:S
        for j = 1:S

            if i < j

                p1 = rand;

                if p1 < C_eff(j)

                    K(i,j) = 1;
                    K(j,i) = 0;

                else

                    K(i,j) = 0;
                    K(j,i) = 0;

                end

            end

        end
    end


    % -----------------------------------------------------
    % Community matrix
    % -----------------------------------------------------

    Moff = zeros(S);

    for i = 1:S
        for j = 1:S

            if i < j && K(i,j)==1 && K(j,i)==0

                Moff(i,j) = -abs(normrnd(0,sigma,1,1));
                Moff(j,i) =  abs(normrnd(0,sigma,1,1));

            end

        end
    end


    % -----------------------------------------------------
    % Eigenvalues
    % -----------------------------------------------------

    lambda_off = eig(Moff);


    % -----------------------------------------------------
    % Different d
    % -----------------------------------------------------

    for k = 1:length(val_d)

        d = val_d(k);

        lambda_M = lambda_off + (1-d);

        rho = max(abs(lambda_M));

        J_case2(k,t) = -log(rho);

    end

end


%% =========================================================
%% 4. Cascade + intervality + broad degree distribution
%% =========================================================

for t = 1:n


    % -----------------------------------------------------
    % Niche values
    % -----------------------------------------------------

    eta = zeros(1,S);

    for i = 1:S

        eta(i) = unifrnd(0,1);

    end

    eta = sort(eta);


    % -----------------------------------------------------
    % Niche radius
    % -----------------------------------------------------

    beta = zeros(1,S);
    r = zeros(1,S);

    for i = 1:S

        beta(i) = betarnd(1,1/C-1);

        r(i) = eta(i)*beta(i);

    end


    % -----------------------------------------------------
    % Number of prey species
    % -----------------------------------------------------

    L = zeros(1,S);

    for i = 1:S

        L(i) = binornd(i-1,r(i));

    end


    % -----------------------------------------------------
    % Starting point of feeding interval
    % -----------------------------------------------------

    s = zeros(1,S);

    for i = 1:S

        s(i) = unidrnd(i-L(i));

    end


    % -----------------------------------------------------
    % Generate adjacency matrix K
    % -----------------------------------------------------

    K = zeros(S);

    for i = 1:S
        for j = 1:S

            if j >= s(i) && j <= s(i)+L(i)-1

                K(i,j) = 1;

            else

                K(i,j) = 0;

            end

        end
    end


    % -----------------------------------------------------
    % Community matrix
    % -----------------------------------------------------

    Moff = zeros(S);

    for i = 1:S
        for j = 1:S

            if K(i,j)==1 && K(j,i)==0

                Moff(i,j) =  abs(normrnd(0,sigma,1,1));
                Moff(j,i) = -abs(normrnd(0,sigma,1,1));

            end

        end
    end


    % -----------------------------------------------------
    % Eigenvalues
    % -----------------------------------------------------

    lambda_off = eig(Moff);


    % -----------------------------------------------------
    % Different d
    % -----------------------------------------------------

    for k = 1:length(val_d)

        d = val_d(k);

        lambda_M = lambda_off + (1-d);

        rho = max(abs(lambda_M));

        J_case3(k,t) = -log(rho);

    end

end


%% =========================================================
%% Mean stability
%% =========================================================

J_casmean = mean(J_cas,2);

J_1mean = mean(J_case1,2);
J_2mean = mean(J_case2,2);
J_3mean = mean(J_case3,2);


%% =========================================================
% Difference relative to the original cascade food web
%
% Delta Stability
% =
% Stability_variant - Stability_cascade
%% =========================================================

Delta = zeros(3,3);

% d = 0.5
Delta(1,1) = J_1mean(1) - J_casmean(1);
Delta(1,2) = J_2mean(1) - J_casmean(1);
Delta(1,3) = J_3mean(1) - J_casmean(1);

% d = 1
Delta(2,1) = J_1mean(2) - J_casmean(2);
Delta(2,2) = J_2mean(2) - J_casmean(2);
Delta(2,3) = J_3mean(2) - J_casmean(2);

% d = 1.5
Delta(3,1) = J_1mean(3) - J_casmean(3);
Delta(3,2) = J_2mean(3) - J_casmean(3);
Delta(3,3) = J_3mean(3) - J_casmean(3);


%% =========================================================
% Same y-axis range for all three panels
%% =========================================================

all_delta = Delta(:);

ymin = min([0;all_delta]);
ymax = max([0;all_delta]);

yrange = ymax-ymin;

if yrange < 1e-10

    yrange = 0.1;

end

ylim_common = [ymin-0.15*yrange, ymax+0.05*yrange];


%% =========================================================
%% Plot
%% =========================================================

figure

set(gcf,...
    'unit','centimeters',...
    'position',[10 5 20 6]);


category_names = {...
    'Intervality',...
    'Broad degree dist.',...
    'Intervality + broad degree dist.'};


%% =========================================================
% d = 0.5
%% =========================================================

subplot(1,3,1)

set(gca,...
    'unit','centimeters',...
    'position',[1.5 1.5 5 3.5],...
    'FontSize',8);

hold on


bar(...
    1:3,...
    Delta(1,:),...
    0.3,...
    'FaceColor',[0 0.4470 0.7410],...
    'EdgeColor','k',...
    'LineWidth',0.6);


yline(0,...
    '-',...
    'Color',[0.3 0.3 0.3],...
    'LineWidth',0.7);


xlim([0.5 3.5]);
ylim(ylim_common);

% 关键：负值在上，正值在下
set(gca,'YDir','reverse');

xticks(1:3);
xticklabels(category_names);

xtickangle(20);

ylabel('$$\Delta$$ Stability',...
    'Interpreter','latex',...
    'FontSize',9);

title('$$d=0.5$$',...
    'Interpreter','latex',...
    'FontSize',10);

box on


%% =========================================================
% d = 1
%% =========================================================

subplot(1,3,2)

set(gca,...
    'unit','centimeters',...
    'position',[7.5 1.5 5 3.5],...
    'FontSize',8);

hold on


bar(...
    1:3,...
    Delta(2,:),...
    0.3,...
    'FaceColor',[0 0.4470 0.7410],...
    'EdgeColor','k',...
    'LineWidth',0.6);


yline(0,...
    '-',...
    'Color',[0.3 0.3 0.3],...
    'LineWidth',0.7);


xlim([0.5 3.5]);
ylim(ylim_common);

% 关键：负值在上，正值在下
set(gca,'YDir','reverse');


xticks(1:3);
xticklabels(category_names);

xtickangle(20);

ylabel('$$\Delta$$ Stability',...
    'Interpreter','latex',...
    'FontSize',9);

title('$$d=1$$',...
    'Interpreter','latex',...
    'FontSize',10);

box on


%% =========================================================
% d = 1.5
%% =========================================================

subplot(1,3,3)

set(gca,...
    'unit','centimeters',...
    'position',[13.5 1.5 5 3.5],...
    'FontSize',8);

hold on


bar(...
    1:3,...
    Delta(3,:),...
    0.3,...
    'FaceColor',[0 0.4470 0.7410],...
    'EdgeColor','k',...
    'LineWidth',0.6);


yline(0,...
    '-',...
    'Color',[0.3 0.3 0.3],...
    'LineWidth',0.7);


xlim([0.5 3.5]);
ylim(ylim_common);

% 关键：负值在上，正值在下
set(gca,'YDir','reverse');


xticks(1:3);
xticklabels(category_names);

xtickangle(20);

ylabel('$$\Delta$$ Stability',...
    'Interpreter','latex',...
    'FontSize',9);

title('$$d=1.5$$',...
    'Interpreter','latex',...
    'FontSize',10);

box on