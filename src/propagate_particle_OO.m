function [out, trajectory] = propagate_particle_OO(obj, in, phi)
    
    % strip the phase space
    x0  = in(1);
    y0  = in(2);
    px0 = in(3);
    py0 = in(4);
    
    % get obj.fields
    
    Ex = obj.E0*obj.fields.Ex;
    Ey = obj.E0*obj.fields.Ey;
    Hz = obj.E0*obj.fields.Hz;

    
    v0 = obj.beta*obj.c0;               % initial speed
    prop_time = obj.num_periods*obj.lambda*obj.beta/v0;           % time for a straight particle to traverse structure
    
    dt = prop_time/1000;            % time step
    TMAX = round(prop_time/dt*2);   % max number of times steps to consider
        
    x = x0;                         % starting x position
    y = y0;                         % starting y position
    px = px0;                       % starting x momentum
    py = py0;                       % starting y momentum
        
    out = [0 0 0 0];                % default return value for 'out' vector (if failed)
    trajectory = zeros(4,TMAX);     % trajectory recording array
    
    % loop through times
    for T = (1:TMAX)
        
        t = (T-1)*dt;
        time_phase = exp(1i*obj.w0*t + 1i*phi);
        
        trajectory(:,T) = [x,y,px,py];
        
        xi = mod(round(x/obj.dl),obj.Nx)+1;
        yi = mod(round(y/obj.dl),obj.Ny)+1;
        
        if (xi <= 0 || yi <= 0 || yi > obj.Ny)
            break;
        end
        
        if (xi > obj.Nx*obj.num_periods)
            out = [x y px py];
            break;            
        end   
        
        Fx = obj.q*real(Ex(xi,yi)*time_phase) + obj.q*py/obj.me/obj.c0*real(Hz(xi,yi)*time_phase);
        Fy = obj.q*real(Ey(xi,yi)*time_phase) + obj.q*px/obj.me/obj.c0*real(Hz(xi,yi)*time_phase);                
        
        px = px + Fx*dt;
        py = py + Fy*dt;
        
        x = x + px*dt/obj.me;
        y = y + py*dt/obj.me; 
                
    end    
    
    trajectory = trajectory(:,1:T);    
end