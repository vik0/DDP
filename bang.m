% returns a vector normlaised to the max possible value
function v_out = bang(v,v_max)
    % nommalize when norm is not near 0
    if norm(v) < 0.9*v_max
        v_out = v;
    else 
        v_out =  v_max*v/norm(v);
    end
end