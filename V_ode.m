% ODE for V function

function dV = V_ode(t,V,xy,g,c_arg)

    % calculate Vx,Vxx
    Vx(1:4,1)    =  V(2:5);
    Vxx(1:4,1:4) = [V(6:9), V(10:13),V(14:17),V(18:21)];
    
    [Qx,Qu,Qv,Qxx,Quu,Qvv,Qux,Qvx,Quv,L,~,~,~] = Q_function(Vx,Vxx,xy,g,c_arg);
    Qvu = Quv; 
    
    Iu = -pinv(Quu - Quv*(Qvv\Qvu))*(Qu - Quv*(Qvv\Qv));
    Iv = -pinv(Qvv - Qvu*(Quu\Quv))*(Qv - Qvu*(Quu\Qu));
    Ku = -pinv(Quu - Quv*(Qvv\Qvu))*(Qux- Quv*(Qvv\Qvx));
    Kv = -pinv(Qvv - Qvu*(Quu\Quv))*(Qvx- Qvu*(Quu\Qux));

    dV(1) = -( L + Iu'*Qu + Iv'*Qv + 0.5*Iu'*Quu*Iu + Iu'*Quv*Iv + 0.5*Iv'*Qvv*Iv);
     
    dV(2:5) = -(Qx + Ku'*Qu+ Kv'*Qv + Qux'*Iu + Qvx'*Iv + Ku'*Quu*Iu +  ...
                    Ku'*Quv*Iv + Kv'*Qvu*Iu + Kv'*Qvv*Iv);
    ddV = -(2*Ku'*Qux + 2*Kv'*Qvx + 2*Kv'*Qvu*Ku + Ku'*Quu*Ku + Kv'*Qvv*Kv + Qxx);
    
    dV(6:21) = ddV(:); % column wise stacking
    
    dV = dV';