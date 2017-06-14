% ODE for initial strategy
% Inputs:
% 

function uv = init_action(t,xy,xe_goal,u_max,v_max)
    % evader
    x(1:2,1) = xy(1:2);
    % pursuer
    y(1:2,1) = xy(3:4);
    
    if length(xe_goal)~=2
        error('dimension of goal should be 2');
    end
    
    g(1:2,1) = xe_goal;
    
    u = bang(g-x,u_max);
    v = bang(x-y,v_max);
    
    uv(1:4,1) = [u;v];
    
end