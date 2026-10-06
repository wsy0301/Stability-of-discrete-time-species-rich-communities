%% =========================================================
% Function 1
%
% 给定椭圆 c,a,b，寻找最大模特征值
%
% lambda(theta)
% =
% c + a*cos(theta) + i*b*sin(theta)
%
% ==========================================================

function [rho_max,lambda_dom,dom_type,theta_dom] = ...
    dominant_ellipse(c,a,b)


    %% -----------------------------------------------------
    % 椭圆退化成一个点
    % ------------------------------------------------------

    if abs(a)<1e-14 && abs(b)<1e-14

        lambda_dom=c;

        rho_max=abs(c);

        dom_type="point";

        theta_dom=NaN;

        return

    end



    %% -----------------------------------------------------
    % 如果 a >= b
    %
    % 最大模必然位于实轴端点
    % ------------------------------------------------------

    if a>=b

        lambda_right=c+a;

        lambda_left=c-a;


        if abs(lambda_right)>=abs(lambda_left)

            lambda_dom=lambda_right;

            theta_dom=0;

        else

            lambda_dom=lambda_left;

            theta_dom=pi;

        end


        rho_max=abs(lambda_dom);

        dom_type="endpoint";

        return

    end



    %% -----------------------------------------------------
    % b > a
    %
    % 此时可能是：
    %
    % 1. 斜向弧上的点
    % 2. 实轴端点
    %
    % 弧上极值满足
    %
    % cos(theta*) = a*c/(b^2-a^2)
    %
    % ------------------------------------------------------

    u_star=a*c/(b^2-a^2);



    if abs(u_star)<=1

        %% -------------------------------------------------
        % 最大模位于斜向弧
        % --------------------------------------------------

        theta_dom=acos(u_star);


        x_dom=c+a*u_star;

        y_dom=b*sqrt(max(0,1-u_star^2));


        % 返回虚部为正的共轭特征值
        lambda_dom=x_dom+1i*y_dom;


        rho_max=abs(lambda_dom);

        dom_type="arc";


    else

        %% -------------------------------------------------
        % 弧上的驻点不存在
        % 最大模由实轴端点决定
        % --------------------------------------------------

        lambda_right=c+a;

        lambda_left=c-a;


        if abs(lambda_right)>=abs(lambda_left)

            lambda_dom=lambda_right;

            theta_dom=0;

        else

            lambda_dom=lambda_left;

            theta_dom=pi;

        end


        rho_max=abs(lambda_dom);

        dom_type="endpoint";

    end

end
