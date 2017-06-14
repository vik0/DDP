function uv = uv_interp(t_point,xy,xe_goal,u_max,v_max) 
    global u_new v_new t_bar;
    uv = [bang(interp1(t_bar,u_new,t_point),u_max)'; bang(interp1(t_bar,v_new,t_point),v_max)']; 
    if sum(isnan(xy))
        t_point,xy,xe_goal,u_max,v_max
        error('nan xy in uv_interp')
    end
end