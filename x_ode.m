% ODE for state variable

function dx = x_ode(t,x,g,u_max,v_max,uv_fun)
    uv = uv_fun(t,x,g,u_max,v_max);
    global u_temp;
    u_temp(end+1,:) = uv;
    % single integrator
    dx(1:4,1) = uv;
    if sum(isnan(dx))
        x
        uv
        uv_fun
        error('error in x_ode')
    end
end