function [out, trajectory] = propagate_particle_OO(obj, in, phi)
    
    % strip the phase space
    z0  = in(1);
    y0  = in(2);
    pz0 = in(3);
    py0 = in(4);
    
    % get obj.fields
    
    Ez = obj.E0*obj.fields.Ez;
    Ey = obj.E0*obj.fields.Ey;
    Hx = obj.E0*obj.fields.Hx;
    
    v0 = obj.c0*obj.beta;
    prop_time = obj.num_periods*obj.lambda*obj.beta/v0;         % time for a straight particle to traverse structure
    
    dt = obj.dl/v0/10;  % time step
    TMAX = round(prop_time/dt*2);   % max number of times steps to consider
        
    z = z0;                         % starting x position
    y = y0;                         % starting y position
    pz = pz0;                       % starting x momentum
    py = py0;                       % starting y momentum
    vz = pz0/obj.me/obj.gamma;
    vy = py0/obj.me/obj.gamma;
    
    out = [0 0 0 0];                % default return value for 'out' vector (if failed)
    trajectory = zeros(4,TMAX);     % trajectory recording array
    
    % loop through times
    for T = (1:TMAX)
        
        t = (T-1)*dt;
        time_phase = exp(-1i*obj.w0*t + 1i*phi);
        
        p_total = sqrt(pz^2 + py^2);
        v_total = sqrt(vz^2 + vy^2);
        
        gamma = p_total/obj.me/v_total; 
        
        trajectory(:,T) = [z,y,pz,py];
        
        zi = mod(round(z/obj.dl),obj.Nz)+1;
        yi = round((y/obj.dl)+obj.ny);
        
        if (zi <= 0 || yi <= 0 || yi > obj.Ny)
            break;
        end
        
        if (z > obj.lambda*obj.beta*obj.num_periods)
            out = [z y pz py];
            break;            
        end           
        
        Fz = obj.q*real(Ez(zi,yi)*time_phase) + obj.q*py/gamma/obj.me/obj.c0*real(Hx(zi,yi)*time_phase);
        Fy = obj.q*real(Ey(zi,yi)*time_phase) + obj.q*pz/gamma/obj.me/obj.c0*real(Hx(zi,yi)*time_phase);                
        
        pz = pz + Fz*dt;
        py = py + Fy*dt;
              
        vz = pz/obj.me/gamma;
        vy = py/obj.me/gamma;
        
        z = z + vz*dt;
        y = y + vy*dt; 
                
    end    
    
    trajectory = trajectory(:,1:T);    
end