close all
clear all
clc

% starting state
x_init = [10;10;10;-10];

global x1_goal;
x1_goal = [0;0];

global u_max v_max;
u_max = 1; v_max = 0.5;

% global varible for storing control actions
global u_temp;

% time horizon
t0 = 0;   tf = 15;

% Calculate the initial trajcetory
tic;    
[t, x_bar] = ode45(@(t,x)x_ode(t,x,x1_goal,u_max,v_max,@init_action),[t0,tf],x_init);
u_bar = u_temp;
u_temp = [];
[t,uniq,~] = unique(t,'stable');
x_bar = x_bar(uniq,1:4);
u_bar = u_bar(uniq,1:4);
toc;

figure;
plot(x_bar(:,1),x_bar(:,2),'r'); hold on;
plot(x_bar(:,3),x_bar(:,4),'b'); 
plot(x1_goal(1),x1_goal(2),'+r');
title('Initial trajectories');
xlabel('x'); ylabel('y');
grid on; axis('equal');
legend('evader','pursuer','goal')

% --- goal has been perturbed ---
x1_goal = [-5;5];
obst = [1,-1];

x_bar_i = x_bar;
u_bar_i = u_bar;
t_i = t;

max_iter =10;
tic;
c_arg = 1; 
for iter = 1:max_iter
    
    % -- terminal value function ---
    [V_f,dV_f,ddV_f]= Val(x_bar_i(end,1:4),x1_goal);
    V_tf = [V_f;dV_f;ddV_f(:)]; % columns of hessian are stacked
    if sum(isnan(V_tf(:)))
        error('nan V_tf')
    end
    clear V_f dV_f ddV_f uniq;

    % ----  propogate V function backward  ------
    V_bar = zeros(length(t_i),length(V_tf));
    V_bar(end,:) = V_tf;
    %     tic;
    for i = length(t_i):-1:2
%         [ti,Vi] = ode45(@(t,V)V_ode(t,V,x_bar_i(i,:),x1_goal,c_arg),[t_i(i),t_i(i-1)],V_tf);
        [ti,Vi] = ode45(@(t,V)V_ode_obst(t,V,x_bar_i(i,:),x1_goal,c_arg,obst),[t_i(i),t_i(i-1)],V_tf);
        V_bar(i-1,:) = Vi(end,:);
        V_tf = Vi(end,:);
    end
    %     toc;

    % ---- computing du and dv 
    dx = []; du = []; dv = []; 
    dx_init = zeros(4,1);
    dx(1,:) = zeros(4,1);
    du(1,:) = zeros(2,1);
    dv(1,:) = zeros(2,1);

    %     tic;
    for i = 1:length(t_i)-1
        Vx(1:4,1) = V_bar(i,2:5)';
        if isnan(Vx)
            error('nan Vx');
        end
        Vxx = reshape(V_bar(i,6:21),[4,4]);
        [Qx,Qu,Qv,Qxx,Quu,Qvv,Qux,Qvx,Quv,L,Fx,Fu,Fv] = Q_function_obst(Vx,Vxx,x_bar_i(i,:),x1_goal,c_arg,obst);
        Qvu = Quv;

        Iu = -pinv(Quu - Quv*(Qvv\Qvu))*(Qu - Quv*(Qvv\Qv));
        Iv = -pinv(Qvv - Qvu*(Quu\Quv))*(Qv - Qvu*(Quu\Qu));
        Ku = -pinv(Quu - Quv*(Qvv\Qvu))*(Qux- Quv*(Qvv\Qvx));
        Kv = -pinv(Qvv - Qvu*(Quu\Quv))*(Qvx- Qvu*(Quu\Qux));
    
        if sum(isnan(Iu(:))) + sum(isnan(Iv(:))) +sum(isnan(Ku(:))) +sum(isnan(Kv(:)))
            error('nan in Iu Iv Ku Kv')
        end
        [~,dxi] = ode45(@(t,dx)dx_ode(t,dx,Fx,Fu,Fv,Iu,Iv,Ku,Kv),[t_i(i),t_i(i+1)],dx_init);
        if sum(isnan(dxi(:)))
            error('nan in dxi');
        end
        dx_init(1:4,1) = dxi(end,:);
        dx(end+1,:) = dx_init;
        du(end+1,:) = Iu + Ku*dx_init;
        dv(end+1,:) = Iv + Kv*dx_init;
    end
    %     toc;
    
    if sum(isnan(du))+ sum(isnan(dv))
        error('nan in du and dv')
    end
    
    gamma = 0.1;
    global t_bar u_new v_new;
    t_bar = (1-gamma)*t_i + gamma*t_i;
    u_new = (1-gamma)*u_bar_i(:,1:2) +gamma* du;
    v_new = (1-gamma)*u_bar_i(:,3:4) +gamma* dv;
    x_bar2 = []; t2 = [];
    % new states from obtained control actions
    [t2, x_bar2] = ode45(@(t,x)x_ode(t,x,x1_goal,u_max,v_max,@uv_interp),[t0,tf],x_init);
    if sum(isnan(x_bar2(:)))
        error('nan in x_bar2');
    end
    u_bar2 = u_temp;
    u_temp = [];
    [t2,uniq,~] = unique(t2,'stable');
    x_bar2 = x_bar2(uniq,1:4);
    u_bar2 = u_bar2(uniq,1:4);
    u_bar_i = u_bar2;
    x_bar_i = x_bar2;
    t_i = t2;
    toc;
    
    
%     figure(2);
%     plot(t,u_bar(:,1),'--g'); hold on;
%     plot(t,u_bar(:,2),'g'); hold on;
%     plot(t,u_bar(:,3),'--c'); hold on;
%     plot(t,u_bar(:,4),'c'); hold on;
%     plot(t,u_new(:,1),'--r'); hold on;
%     plot(t,u_new(:,2),'r'); hold on;
%     plot(t,v_new(:,1),'--b');
%     plot(t,v_new(:,2),'b');
%     legend('u_bar_x','u_bar_y','v_bar_x','v_bar_y','u_new_x','u_new_y','v_new_x','v_new_y'); 
%     grid on; axis('equal');


    figure(11);
    plot(x_bar(:,1),x_bar(:,2),'--r'); hold on;
    plot(x_bar2(:,1),x_bar2(:,2),'r'); hold on;
    plot(x_bar(:,3),x_bar(:,4),'--b');
    plot(x_bar2(:,3),x_bar2(:,4),'b');
    grid on; axis ('equal'); hold off;
    legend('e_{old}','e_{new}','p_{old}','p_{new}');
    pause(0.1);
end

toc;
figure;
plot(x_bar(:,1),x_bar(:,2),'--r'); hold on;
plot(x_bar2(:,1),x_bar2(:,2),'r'); hold on;
plot(x_bar(:,3),x_bar(:,4),'--b');
plot(x_bar2(:,3),x_bar2(:,4),'b');
grid on; axis ('equal');
legend('e_{old}','e_{new}','p_{old}','p_{new}');