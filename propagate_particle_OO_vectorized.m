function [outs] = propagate_particle_OO_vectorized(obj, ins)
    
    [~,Ne] = size(ins);
    outs = zeros(4,Ne);   
    
    % get obj.fields
    
    Ez = obj.E0*obj.fields.Ez;
    Ey = obj.E0*obj.fields.Ey;
    Hx = obj.E0*obj.fields.Hx;
    
    v0 = obj.c0*obj.beta;
    prop_time = obj.num_periods*obj.lambda*obj.beta/v0;         % time for a straight particle to traverse structure
    
    dt = obj.dl/v0/10;  % time step
    TMAX = round(prop_time/dt*2);   % max number of times steps to consider
    
    curr_vec = ins(1:4,:);
    phi_vec = ins(5,:);
    vel_vec = ins(3:4,:)/obj.me/obj.gamma;
    
    still_running_vec = ones(1,Ne);
    
    % loop through times
    for T = (1:TMAX)
        
        t = (T-1)*dt;
        time_phase = exp(-1i*obj.w0*t);
        time_phase_vec = time_phase*exp(1i*phi_vec);
        
        p_total_vec = sqrt(curr_vec(3,:).^2 + curr_vec(4,:).^2);
        v_total_vec = sqrt(vel_vec(1,:).^2 + vel_vec(2,:).^2);
        
        gamma_vec = p_total_vec./v_total_vec./obj.me; 
                
        zi_vec = mod(round(curr_vec(1,:)./obj.dl),obj.Nz)+1;
        yi_vec = round((curr_vec(1,:)./obj.dl)+obj.ny);
                       

        still_running_vec(zi_vec<=0) = 0;
        still_running_vec(yi_vec<=0) = 0;
        still_running_vec(yi_vec>obj.Ny) = 0;
        still_running_vec(zi_vec>obj.Nz*obj.num_periods) = 0;

       % outs(zi_vec>obj.Nz*obj.num_periods) = curr_vec(:,zi_vec>obj.Nz*obj.num_periods);
        outs = curr_vec;
        
        zi_vec(still_running_vec==0) = 1;
        yi_vec(still_running_vec==0) = 1;
                
        pos_index_vec = sub2ind([obj.Nz*obj.num_periods,obj.Ny],zi_vec,yi_vec);
        pos_index_vec = mod(pos_index_vec,obj.Nz)+ones(1,Ne);
        
        Ez_vec = Ez(pos_index_vec);
        Ey_vec = Ey(pos_index_vec);
        Hx_vec = Hx(pos_index_vec);
        
        complex_Fz_vec = obj.q*( Ez_vec + curr_vec(4,:)./gamma_vec./obj.me./obj.c0.*Hx_vec);
        complex_Fy_vec = obj.q*( Ey_vec + curr_vec(3,:)./gamma_vec./obj.me./obj.c0.*Hx_vec);
        
        Fz_vec = real(complex_Fz_vec.*time_phase_vec).*still_running_vec;
        Fy_vec = real(complex_Fy_vec.*time_phase_vec).*still_running_vec;

        curr_vec(3,:) = curr_vec(3,:) + Fz_vec*dt;
        curr_vec(4,:) = curr_vec(4,:) + Fy_vec*dt;
              
        vel_vec(1,:) = curr_vec(3,:)./obj.me./gamma_vec;
        vel_vec(2,:) = curr_vec(4,:)./obj.me./gamma_vec;
        
        curr_vec(1,:) = curr_vec(1,:) + vel_vec(1,:)*dt;
        curr_vec(2,:) = curr_vec(2,:) + vel_vec(2,:)*dt;
           
        if ~sum(still_running_vec)
            break
        end
    end    
end