% Function to give Q and F derivatives

function [Qx,Qu,Qv,Qxx,Quu,Qvv,Qux,Qvx,Quv,L,Fx,Fu,Fv] = Q_function(Vx,Vxx,xy,g,funct)

    x(1:2,1) = xy(1:2);
    y(1:2,1) = xy(3:4);
    
    g(1:2,1) = g;
    
    global u_max;
    global v_max;
    
    u = bang(g-x,u_max);
    v = bang(x-y,v_max);
   
    f1 = norm(x-g)^2;
    df1(1:4,1) = 0;
    df1(1:2,1) = 2*(x-g);
    
    ddf1(1:4,1:4) = 0;
    ddf1(1:2,1:2) =  2*eye(2);
    if funct == 1
        f2 = -norm(x-y)^2;
        df2(1:2,1) = -2*(x-y);
        df2(3:4,1) = -2*(y-x);

        ddf2(1:4,1:4) = 0;
        ddf2(1:2,1:2) =  -2*eye(2);
        ddf2(3:4,3:4) =  -2*eye(2);
        ddf2(1:2,3:4) =   2*eye(2);
        ddf2(3:4,1:2) =   2*eye(2);
        
        k2 = 1;
    else
        f2 = 1/(norm(x-y)^2);
        
        df2(1:2,1) = -2*(x-y)*f2^2;
        df2(3:4,1) = -2*(y-x)*f2^2;
        
        ddf2(1:4,1:4) = 0;
        t1 = eye(2)*(-2*f2^2);
        t2 = 8*(x-y)*(x-y)'*f2^3;
        
        %         size(t2)
        
        ddf2(1:2,1:2) = t1 + t2;
        ddf2(3:4,3:4) = t1 + t2;
        ddf2(1:2,3:4) =-t1 + t2;
        ddf2(3:4,1:2) =-t1 + t2;
        
        k2 = 1;
    end
    k = 100;
    
    f = k*(f1)+k2*f2;
    df = k*(df1)+k2*df2;
    ddf = k*(ddf1)+k2*ddf2;
    
    ru = 10; 
    rv = 10;
    L   =  f+ ru*norm(u)^2 - rv*norm(v)^2;
    Lx  =  df;
    Lxx =  ddf;
    Lu  =  2*ru*u ; Lv = -2*rv* v;
    Luu =  2*ru*eye(2);
    Lvv = -2*rv*eye(2);
    Luv = zeros(2);
    Lux = zeros(2,4);
    Lvx = zeros(2,4);              
     
    Fx(1:4,1:4) = zeros(4) + rand*0.001*eye(4);
    Fu(1:4,1:2) = [1,0;0,1;0,0;0,0];
    Fv(1:4,1:2) = [0,0;0,0;1,0;0,1];
    
    Qx = Fx'*Vx + Lx;
    Qu = Fu'*Vx + Lu;
    Qv = Fv'*Vx + Lv;
    Qxx = Lxx + 2*Vxx*Fx + rand*0.001*eye(size(Lxx));
    Quu = Luu + rand*0.001*eye(size(Luu));
    Qvv = Lvv + rand*0.001*eye(size(Lvv));
    Qux = Fu'*Vxx + Lux;
    Qvx = Fv'*Vxx + Lvx;
    Quv = Luv;
    
end