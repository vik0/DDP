% ODE for delta_x
function ddx = dx_ode(t,dx,Fx,Fu,Fv,Iu,Iv,Ku,Kv) 
    Iu(abs(Iu)>1e4) = 1e4*sign(Iu(abs(Iu)>1e4));
    Iv(abs(Iv)>1e4) = 1e4*sign(Iv(abs(Iv)>1e4));
    Ku(abs(Ku)>1e4) = 1e4*sign(Ku(abs(Ku)>1e4));
    Kv(abs(Kv)>1e4) = 1e4*sign(Kv(abs(Kv)>1e4));
    
    ddx =  (Fu*Iu + Fv*Iv + (Fx + Fu*Ku + Fv*Kv)*dx);
    if sum(isnan(ddx))
        
        dx,Fx,Fu,Fv,Iu,Iv,Ku,Kv, ddx,dx
        error('nan in ddx');
    end
end