% Function returns value function, its gradient, hessian

% TODO: Implement more complex functions

function [f,df,ddf] = Val(xy,goal)
    
    x(1:2,1) = xy(1:2);
    y(1:2,1) = xy(3:4);
    g(1:2,1) = goal;
    
    k = 100;
    
    f = norm(x-g)^2;
    df(1:4,1) = 0;
    df(1:2,1) = 2*(x-g);
    
    
    ddf(1:4,1:4) = 0;
    ddf(1:2,1:2) = 2*eye(2);
    
    f2 = norm(x-y)^2;
    df2(1:2,1) = 2*(x-y);
    df2(3:4,1) = 2*(y-x);
    
    ddf2(1:4,1:4) = 0;
    ddf2(1:2,1:2) =  2*eye(2);
    ddf2(3:4,3:4) =  2*eye(2);
    ddf2(1:2,3:4) = -2*eye(2);
    ddf2(3:4,1:2) = -2*eye(2);
    
    f = k*(f)-f2;
    df = k*(df) -df2;
    ddf = k*(ddf) -ddf2;
    
end