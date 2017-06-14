function [f2,df2,ddf2] = cost_der(x0,y0,c_arg,k)
    f2 = 0;
    df2 = zeros(4,1);
    ddf2 = zeros(4);
    
    x(1:2,1) =x0;
    y(1:2,1) =y0;
    
    if nargin ==3
        k = 1;
    end
    
    if c_arg == 1    
        % -ve distance from opponent
        f2 = -k*norm(x-y)^2;
        df2(1:2,1) = -2*k*(x-y);
        df2(3:4,1) = -2*k*(y-x);

        ddf2(1:4,1:4) = 0;
        ddf2(1:2,1:2) =  -2*k*eye(2);
        ddf2(3:4,3:4) =  -2*k*eye(2);
        ddf2(1:2,3:4) =   2*k*eye(2);
        ddf2(3:4,1:2) =   2*k*eye(2);
        
    elseif c_arg == 2
        
        % inverse of distance from opponent 
        f2 = k/(norm(x-y)^2+1e-2);
        df2(1:2,1) = -k*2*(x-y)*f2^2;
        df2(3:4,1) = -k*2*(y-x)*f2^2;
        
        ddf2(1:4,1:4) = 0;
        t1 = k*eye(2)*(-2*f2^2);
        t2 = k*8*(x-y)*(x-y)'*f2^3;
        
        %         size(t2)
        
        ddf2(1:2,1:2) = t1 + t2;
        ddf2(3:4,3:4) = t1 + t2;
        ddf2(1:2,3:4) =-t1 + t2;
        ddf2(3:4,1:2) =-t1 + t2;
       
    elseif c_arg == 3
        if nargin ==3
            k = 0.01;
        end
        f2 = exp(-k*norm(x-y)^2);
        df2(1:2,1) = -k*2*(x-y)*f2;
        df2(3:4,1) = -k*2*(y-x)*f2;
        
        ddf2(1:4,1:4) = 0;
        t1 = -2*k*eye(2)*f2;
        t2 = 4*k*k*(x-y)*(x-y)'*f2;
        %         size(t2)
        
        ddf2(1:2,1:2) =  t1 + t2;
        ddf2(1:2,3:4) = -t1 - t2;
        ddf2(3:4,1:2) = -t1 - t2;
        ddf2(3:4,3:4) =  t1 + t2;    
    end
    
    if f2 == 0
        disp('something is wrong in computing cost');
        x 
        y 
        disp(c_arg)
    end
end