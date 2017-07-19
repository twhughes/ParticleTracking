function [out, trajectory] = propagate_particle(in, fields, num_periods, beta, n_per_lam, E0, phi)
    
    lambda = 2e-6;
    c0 = 3e8;
    q = 1.6e-14;
    me = 9e-31;
    % strip the phase space
    x0  = in(1);
    y0  = in(2);
    px0 = in(3);
    py0 = in(4);
    
    % get fields
    Ex = [];
    Ey = [];
    Hz = [];
    for i = (1:num_periods)
        Ex = [Ex; E0*fields.Ex];
        Ey = [Ey; E0*fields.Ey];
        Hz = [Hz; E0*fields.Hz];        
    end

    [Nx,Ny] = size(Ex);
    dl = lambda/n_per_lam;
    Lx = lambda*num_periods;
    c0 = 3e8;
    v0 = beta*c0;
    t0 = Lx/v0;
    
    dt = t0/100;    
    TMAX = round(t0/dt*20);
        
    x = x0;
    y = y0;
    px = px0;
    py = py0;
    
    out = [0 0 0 0];
    trajectory = zeros(4,TMAX);
    
    for T = (1:TMAX)
               
        time_phase = exp(1i*2*pi*c0/lambda*T*dt + 1i*2*pi*phi);
        
        trajectory(:,T) = [x,y,px,py];
        
        xi = round(x/dl);
        yi = round(y/dl);
        
        if (xi <= 0 || yi <= 0 || yi > Ny)
            break;
        end
        
        if (xi > Nx)
            out = [x y px py];
            break;            
        end   
        
        Fx = q*real(Ex(xi,yi)*time_phase) + q*py/me/c0*real(Hz(xi,yi)*time_phase);
        Fy = q*real(Ey(xi,yi)*time_phase) + q*px/me/c0*real(Hz(xi,yi)*time_phase);
        
        px = px + Fx*dt;
        py = py + Fy*dt;
        
        x = x + px*dt/me;
        y = y + py*dt/me; 
                
    end    
    
    trajectory = trajectory(:,1:T);
    
end