
% This file can be used to generate Fig. 4c 



clear all
clc
close all


%% =========================================================
% Multi-species Ricker model
%
% Only competitive community
%
% Model:
%
% X_i(t+1)
% =
% X_i(t)*exp[
%     r_i
%     - s_i X_i(t)
%     + sum_{j~=i} A_ij X_j(t)
% ]
%
%
% Equilibrium:
%
%       X_i^* = 1
%
%
% At X*=1:
%
%       F_ii = -s_i
%
%       F_ij = A_ij
%
%
% Define:
%
%       F_ii = -d
%
%
% Therefore:
%
%       s_i = d
%
%
% To guarantee X*=1:
%
%       r_i
%       =
%       d - sum_{j~=i} F_ij
%
%
% Interspecific interaction:
%
%       A_ij = F_ij
%
%
% Community matrix:
%
%       M = I + F
%
%
% Stability:
%
%       rho = max(abs(eig(M)))
%
%       alpha = -log(rho)
%
%% =========================================================



%% =========================================================
% Basic parameters
%% =========================================================

% Number of species
S = 100;


% Connectance
C = 0.10;


% Standard deviation of interaction strengths
sigma = 0.01;


% Four self-regulation strengths
d_values = ...
    [0.6 1 1.2 1.6];


n_d = ...
    length(d_values);


% Random seed
rng(1);



%% =========================================================
% Generate competitive Foff
%
% For each species pair:
%
%       probability = C
%
% If connected:
%
%       F_ij < 0
%
%       F_ji < 0
%
%
% Importantly:
%
%       F_ij
%
% and
%
%       F_ji
%
% are independently generated.
%% =========================================================

Foff = ...
    zeros(S,S);


for i = 2:S

    for j = 1:i-1


        p1 = ...
            rand;


        if p1 <= C


            Foff(i,j) = ...
                -abs(normrnd(0,sigma));


            Foff(j,i) = ...
                -abs(normrnd(0,sigma));


        end


    end

end


% Diagonal = 0
Foff(1:S+1:end) = 0;



%% =========================================================
% Realized connectance
%% =========================================================

C_real = ...
    nnz(Foff)/(S*(S-1));


fprintf('\n')

disp('============================================')
disp('Competitive Ricker community')
disp('============================================')


fprintf(...
    'S = %d\n',...
    S);


fprintf(...
    'Target C = %.4f\n',...
    C);


fprintf(...
    'Realized C = %.4f\n',...
    C_real);


fprintf(...
    'sigma = %.4f\n',...
    sigma);



%% =========================================================
% Row sum of off-diagonal feedbacks
%
%       sumF_i
%       =
%       sum_{j~=i} F_ij
%% =========================================================

sumF = ...
    sum(Foff,2);



%% =========================================================
% Interspecific interaction matrix
%
% For Ricker:
%
%       A_ij = F_ij
%% =========================================================

A = ...
    Foff;



%% =========================================================
% Initialize variables
%% =========================================================

F = ...
    cell(1,n_d);


M = ...
    cell(1,n_d);


r = ...
    cell(1,n_d);


s = ...
    cell(1,n_d);


eig_value = ...
    cell(1,n_d);


rho = ...
    zeros(1,n_d);


alpha = ...
    zeros(1,n_d);


equilibrium_error = ...
    zeros(1,n_d);



%% =========================================================
% Calculate Ricker parameters and local stability
%% =========================================================

for kk = 1:n_d


    d = ...
        d_values(kk);



    %% =====================================================
    % Full feedback matrix F
    %
    % Off-diagonal:
    %
    %       F_ij = Foff_ij
    %
    % Diagonal:
    %
    %       F_ii = -d
    %% =====================================================

    F{kk} = ...
        Foff;


    F{kk}(1:S+1:end) = ...
        -d;



    %% =====================================================
    % Community matrix
    %
    %       M = I + F
    %% =====================================================

    M{kk} = ...
        eye(S) + F{kk};



    %% =====================================================
    % Self-regulation
    %
    %       s_i = d
    %% =====================================================

    s{kk} = ...
        d*ones(S,1);



    %% =====================================================
    % Intrinsic growth rate
    %
    % At X*=1:
    %
    %       0
    %       =
    %       r_i
    %       - d
    %       + sum_j F_ij
    %
    % Therefore:
    %
    %       r_i
    %       =
    %       d - sum_j F_ij
    %% =====================================================

    r{kk} = ...
        d - sumF;



    %% =====================================================
    % Eigenvalues
    %% =====================================================

    eig_value{kk} = ...
        eig(M{kk});



    %% =====================================================
    % Spectral radius
    %% =====================================================

    rho(kk) = ...
        max(abs(eig_value{kk}));



    %% =====================================================
    % Stability
    %
    %       alpha = -log(rho)
    %% =====================================================

    alpha(kk) = ...
        -log(rho(kk));



    %% =====================================================
    % Check equilibrium X*=1
    %% =====================================================

    Xstar = ...
        ones(S,1);


    Xcheck = ...
        Xstar .* ...
        exp(...
        r{kk} ...
        - s{kk}.*Xstar ...
        + A*Xstar);


    equilibrium_error(kk) = ...
        max(abs(Xcheck-Xstar));


end



%% =========================================================
% Display local stability
%% =========================================================

disp(' ')

disp('============================================')
disp('LOCAL STABILITY')
disp('============================================')


for kk = 1:n_d


    fprintf(...
        ['d = %.3f   ',...
         'rho = %.6f   ',...
         'alpha = %.6f   ',...
         'equilibrium error = %.2e\n'],...
         d_values(kk),...
         rho(kk),...
         alpha(kk),...
         equilibrium_error(kk));


end



%% =========================================================
% Simulation settings
%% =========================================================

% Maximum simulation time
Tmax = 200;


% Number of time steps displayed
Tshow = 40;


% Extinction threshold
extinction_threshold = ...
    0.0001;




%% =========================================================
% Initial condition
%
% Fixed total perturbation strength
% Random perturbation direction
%% =========================================================

Xstar = ...
    ones(S,1);


% Prescribed total perturbation strength
D0_fixed = ...
    0.5;


% Generate a random direction
direction = ...
    randn(S,1);


% Normalize the direction
direction = ...
    direction / norm(direction,2);


% Initial abundance
X0 = ...
    Xstar + D0_fixed*direction;


% Check the actual perturbation strength
D0 = ...
    norm(X0-Xstar,2);


fprintf(...
    'Prescribed perturbation strength = %.8f\n',...
    D0_fixed);


fprintf(...
    'Actual perturbation strength     = %.8f\n',...
    D0);



%% =========================================================
% Initial perturbation strength
%
% Euclidean distance:
%
%       D0
%       =
%       ||X(0)-X*||_2
%% =========================================================

D0 = ...
    norm(X0-Xstar,2);



%% =========================================================
% Recovery criterion
%
% The community is considered recovered when:
%
%       D(t)
%       <=
%       0.005*D0
%
%
% That is:
%
%       D(t)/D0 <= 0.005
%
%
% The remaining perturbation is only 0.5%
% of the initial perturbation.
%% =========================================================

recovery_fraction = ...
    0.005;


recovery_threshold = ...
    recovery_fraction*D0;



fprintf('\n')

disp('============================================')
disp('RECOVERY CRITERION')
disp('============================================')


fprintf(...
    'Initial distance D0 = %.8f\n',...
    D0);


fprintf(...
    'Recovery threshold = %.8f\n',...
    recovery_threshold);


fprintf(...
    'Threshold = %.2f%% of initial perturbation\n',...
    100*recovery_fraction);



%% =========================================================
% Initialize trajectories
%% =========================================================

X = ...
    cell(1,n_d);


distance = ...
    cell(1,n_d);


recovery_time = ...
    NaN(1,n_d);



%% =========================================================
% Simulate Ricker dynamics
%% =========================================================

for kk = 1:n_d


    %% -----------------------------------------------------
    % Initialize abundance
    %% ------------------------------------------------------

    X{kk} = ...
        zeros(S,Tmax+1);


    X{kk}(:,1) = ...
        X0;



    %% -----------------------------------------------------
    % Initialize Euclidean distance
    %% ------------------------------------------------------

    distance{kk} = ...
        NaN(1,Tmax+1);


    distance{kk}(1) = ...
        D0;



    %% =====================================================
    % Time evolution
    %% =====================================================

    for tt = 1:Tmax


        Xold = ...
            X{kk}(:,tt);



        %% -------------------------------------------------
        % Ricker exponent
        %% --------------------------------------------------

        exponent = ...
            r{kk} ...
            - s{kk}.*Xold ...
            + A*Xold;



        %% -------------------------------------------------
        % Ricker update
        %% --------------------------------------------------

        Xnew = ...
            Xold .* exp(exponent);



        %% -------------------------------------------------
        % Avoid numerical overflow
        %% --------------------------------------------------

        if any(~isfinite(Xnew))


            warning(...
                ['Non-finite abundance at ',...
                 'd = %.3f, t = %d'],...
                 d_values(kk),...
                 tt);


            X{kk}(:,tt+1:end) = ...
                NaN;


            break


        end



        %% -------------------------------------------------
        % Extinction rule
        %% --------------------------------------------------

        Xnew(...
            Xnew < extinction_threshold) = 0;



        %% -------------------------------------------------
        % Save abundance
        %% --------------------------------------------------

        X{kk}(:,tt+1) = ...
            Xnew;



        %% -------------------------------------------------
        % Euclidean distance from equilibrium
        %
        %       D(t)
        %       =
        %       ||X(t)-X*||_2
        %% --------------------------------------------------

        distance{kk}(tt+1) = ...
            norm(Xnew-Xstar,2);


    end



    %% =====================================================
    % Recovery time
    %
    % First time satisfying:
    %
    %       D(t) <= 0.005 D0
    %% =====================================================

    index_recovery = ...
        find(...
        distance{kk} <= recovery_threshold,...
        1,...
        'first');



    if ~isempty(index_recovery)


        recovery_time(kk) = ...
            index_recovery-1;


    end


end



%% =========================================================
% Display recovery time
%% =========================================================

disp(' ')

disp('============================================')
disp('RECOVERY TIME')
disp('============================================')


for kk = 1:n_d


    if isnan(recovery_time(kk))


        fprintf(...
            ['d = %.3f   ',...
             'No recovery within %d steps   ',...
             'rho = %.6f   ',...
             'alpha = %.6f\n'],...
             d_values(kk),...
             Tmax,...
             rho(kk),...
             alpha(kk));


    else


        fprintf(...
            ['d = %.3f   ',...
             'T_R = %.0f   ',...
             'rho = %.6f   ',...
             'alpha = %.6f\n'],...
             d_values(kk),...
             recovery_time(kk),...
             rho(kk),...
             alpha(kk));


    end


end



%% =========================================================
% Time vector
%% =========================================================

time_show = ...
    0:Tshow;



%% =========================================================
% Determine common y-axis limits
%
% The first four subplots use the same y limits.
%% =========================================================

all_abundance = [];


for kk = 1:n_d


    Xtemp = ...
        X{kk}(:,1:Tshow+1);


    value_temp = ...
        Xtemp(isfinite(Xtemp));


    all_abundance = ...
        [all_abundance; value_temp];


end



if isempty(all_abundance)


    abundance_lower = 0.5;

    abundance_upper = 1.5;


else


    abundance_lower = ...
        min(all_abundance);


    abundance_upper = ...
        max(all_abundance);


    abundance_range = ...
        abundance_upper ...
        - abundance_lower;



    if abundance_range < 0.05


        abundance_range = ...
            0.05;


    end



    abundance_lower = ...
        abundance_lower ...
        - 0.08*abundance_range;


    abundance_upper = ...
        abundance_upper ...
        + 0.08*abundance_range;


    abundance_lower = ...
        max(0,abundance_lower);


    % Prevent extreme unstable trajectories
    % from making all stable panels unreadable.
    abundance_upper = ...
        min(abundance_upper,3);


end



%% =========================================================
%% ONE FIGURE
%
% First row:
%
%       4 abundance panels
%
% Each panel:
%
%       4.5 cm x 3 cm
%
%
% Second row:
%
%       recovery time + stability
%
% Also:
%
%       4.5 cm x 3 cm
%% =========================================================

figure


set(...
    gcf,...
    'unit','centimeters',...
    'position',[3 2 22 10.5]);



%% =========================================================
%% subplot 1
%
% d = d_values(1)
%% =========================================================

kk = 1;


subplot(2,4,1)


set(...
    gca,...
    'unit','centimeters',...
    'position',[1.2 6.2 4.5 3]);


hold on


Xshow = ...
    X{kk}(:,1:Tshow+1);



%% ---------------------------------------------------------
% All species
%% ---------------------------------------------------------

plot(...
    time_show,...
    Xshow',...
    'LineWidth',0.45);



%% ---------------------------------------------------------
% Equilibrium abundance
%% ---------------------------------------------------------

yline(...
    1,...
    ':',...
    'LineWidth',0.8);



xlabel(...
    'Time',...
    'FontSize',9);


ylabel(...
    'Abundance',...
    'FontSize',9);



xlim([0 Tshow])


ylim(...
    [abundance_lower abundance_upper])


set(...
    gca,...
    'FontSize',8);



title(...
    sprintf(...
    '$d=%.1f$',...
    d_values(kk)),...
    'Interpreter','latex',...
    'FontSize',9);



if isnan(recovery_time(kk))


    txt = ...
        sprintf(...
        '$T_R=\\mathrm{NR},\\;\\alpha=%.3f$',...
        alpha(kk));


else


    txt = ...
        sprintf(...
        '$T_R=%.0f,\\;\\alpha=%.3f$',...
        recovery_time(kk),...
        alpha(kk));


end



text(...
    0.96,...
    0.92,...
    txt,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'VerticalAlignment','top',...
    'Interpreter','latex',...
    'FontSize',7);


box on



%% =========================================================
%% subplot 2
%
% d = d_values(2)
%% =========================================================

kk = 2;


subplot(2,4,2)


set(...
    gca,...
    'unit','centimeters',...
    'position',[6.4 6.2 4.5 3]);


hold on


Xshow = ...
    X{kk}(:,1:Tshow+1);



plot(...
    time_show,...
    Xshow',...
    'LineWidth',0.45);



yline(...
    1,...
    ':',...
    'LineWidth',0.8);



xlabel(...
    'Time',...
    'FontSize',9);



xlim([0 Tshow])


ylim(...
    [abundance_lower abundance_upper])


set(...
    gca,...
    'FontSize',8);



title(...
    sprintf(...
    '$d=%.1f$',...
    d_values(kk)),...
    'Interpreter','latex',...
    'FontSize',9);



if isnan(recovery_time(kk))


    txt = ...
        sprintf(...
        '$T_R=\\mathrm{NR},\\;\\alpha=%.3f$',...
        alpha(kk));


else


    txt = ...
        sprintf(...
        '$T_R=%.0f,\\;\\alpha=%.3f$',...
        recovery_time(kk),...
        alpha(kk));


end



text(...
    0.96,...
    0.92,...
    txt,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'VerticalAlignment','top',...
    'Interpreter','latex',...
    'FontSize',7);


box on



%% =========================================================
%% subplot 3
%
% d = d_values(3)
%% =========================================================

kk = 3;


subplot(2,4,3)


set(...
    gca,...
    'unit','centimeters',...
    'position',[11.6 6.2 4.5 3]);


hold on


Xshow = ...
    X{kk}(:,1:Tshow+1);



plot(...
    time_show,...
    Xshow',...
    'LineWidth',0.45);



yline(...
    1,...
    ':',...
    'LineWidth',0.8);



xlabel(...
    'Time',...
    'FontSize',9);



xlim([0 Tshow])


ylim(...
    [abundance_lower abundance_upper])


set(...
    gca,...
    'FontSize',8);



title(...
    sprintf(...
    '$d=%.1f$',...
    d_values(kk)),...
    'Interpreter','latex',...
    'FontSize',9);



if isnan(recovery_time(kk))


    txt = ...
        sprintf(...
        '$T_R=\\mathrm{NR},\\;\\alpha=%.3f$',...
        alpha(kk));


else


    txt = ...
        sprintf(...
        '$T_R=%.0f,\\;\\alpha=%.3f$',...
        recovery_time(kk),...
        alpha(kk));


end



text(...
    0.96,...
    0.92,...
    txt,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'VerticalAlignment','top',...
    'Interpreter','latex',...
    'FontSize',7);


box on



%% =========================================================
%% subplot 4
%
% d = d_values(4)
%% =========================================================

kk = 4;


subplot(2,4,4)


set(...
    gca,...
    'unit','centimeters',...
    'position',[16.8 6.2 4.5 3]);


hold on


Xshow = ...
    X{kk}(:,1:Tshow+1);



plot(...
    time_show,...
    Xshow',...
    'LineWidth',0.45);



yline(...
    1,...
    ':',...
    'LineWidth',0.8);



xlabel(...
    'Time',...
    'FontSize',9);



xlim([0 Tshow])


ylim(...
    [abundance_lower abundance_upper])


set(...
    gca,...
    'FontSize',8);



title(...
    sprintf(...
    '$d=%.1f$',...
    d_values(kk)),...
    'Interpreter','latex',...
    'FontSize',9);



if isnan(recovery_time(kk))


    txt = ...
        sprintf(...
        '$T_R=\\mathrm{NR},\\;\\alpha=%.3f$',...
        alpha(kk));


else


    txt = ...
        sprintf(...
        '$T_R=%.0f,\\;\\alpha=%.3f$',...
        recovery_time(kk),...
        alpha(kk));


end



text(...
    0.96,...
    0.92,...
    txt,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'VerticalAlignment','top',...
    'Interpreter','latex',...
    'FontSize',7);


box on



%% =========================================================
%% subplot 5
%
% Recovery time and stability
%% =========================================================

subplot(2,4,5)


set(...
    gca,...
    'unit','centimeters',...
    'position',[8.75 1.2 4.5 3]);


hold on



% =========================================================
% x positions
% =========================================================

x_position = ...
    1:4;



% =========================================================
% LEFT AXIS
%
% Recovery time
% =========================================================

yyaxis left


b1 = ...
    bar(...
    x_position-0.14,...
    recovery_time,...
    0.25);



ylabel(...
    'Recovery time',...
    'FontSize',9);



% =========================================================
% Left y-axis limit
% =========================================================

finite_recovery = ...
    recovery_time(...
    ~isnan(recovery_time));



if isempty(finite_recovery)


    recovery_ymax = 10;


else


    recovery_ymax = ...
        max(finite_recovery);


    recovery_ymax = ...
        max(recovery_ymax,1);


end



ylim(...
    [0 1.15*recovery_ymax]);



% =========================================================
% Mark no-recovery cases
% =========================================================

for kk = 1:n_d


    if isnan(recovery_time(kk))


        text(...
            x_position(kk)-0.14,...
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
% =========================================================

yyaxis right


b2 = ...
    bar(...
    x_position+0.14,...
    alpha,...
    0.25);



ylabel(...
    'Stability',...
    'FontSize',9);



% =========================================================
% Stability boundary
% =========================================================

yline(...
    0,...
    ':',...
    'LineWidth',0.5);

ylim([0 3.1])



% =========================================================
% x axis
% =========================================================

d_labels = ...
    {sprintf('%.1f',d_values(1)),...
     sprintf('%.1f',d_values(2)),...
     sprintf('%.1f',d_values(3)),...
     sprintf('%.1f',d_values(4))};



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



% =========================================================
% Legend
% =========================================================

legend(...
    [b1 b2],...
    {'Recovery time','Stability'},...
    'Location','best',...
    'Box','off',...
    'FontSize',6);


title('Multispecies Ricker model','FontSize',9); 


% =========================================================
% Force x-axis to black
%% =========================================================

ax = ...
    gca;


ax.XAxis.Color = ...
    [0 0 0];


ax.XAxis.Label.Color = ...
    [0 0 0];


box on