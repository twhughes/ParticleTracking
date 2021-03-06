classdef simulation < dynamicprops
    % Class used for simulating particle tracking for dielectric laser accelerators.      
    % 
    %     In this documentation, the 'obj' in 'obj.x' refers to the simulation object that is created by obj = simulation.    
    %     The 'x' refers to either the parameter (below) or the method (further below)
    
    properties (Constant)        
        % constants        
        c0 = 299792458;                     % speed of light (m/s)
        e0 = 8.854187817e-12;               % vacuum permittivity (F/m)
        Z0 = 376.730313416;                 % impedance of free space (Ohms)
        q  = 1.6021766208e-19;              % fundamental charge (C)
        me = 9.10938356e-31;                % mass of electron (kg)
    end
    
    properties     
        % Input pulse default parameters
        lambda = 2e-6;                      % wavelength (m)
        w0;                                 % angular frequency (rad/sec)
        E0 = 1e9;                           % input field strength (V/m)
        % Geometry parameters
        beta = 0.5;                         % electron speed / speed of light
        gamma;                              % relativistic gamma factor of FDFD test particle
        num_periods = 20;                   % number of DLA periods to simulate
        r_pillars = 635.26e-9;              % pillar radius (m)
        gap = 400e-9;                       % gap spacing (m)
        L;                                  % length of structure along propagation
        W;                                  % width of structure perpendicular to propagation
        Nz;                                 % number of x grid points
        Ny;                                 % number of y grid points
        nz;                                 % number of x grid points/2
        ny;                                 % number of y grid points/2
        % Field simulation parameters
        n_per_lam = 200;                    % number of grid points per wavelenth
        dual_drive = 1;                     % dual sided drive (make 0 if single side, 1 if dual side)
        npml = 50;                          % number of PML points
        eps = 3.45^2;                       % relative permittivity of pillars        
        dl;
        % DLA parameters
        G;                                  % gradient of one section (V/m)   
        phi0;                               % phase of light for maximum acceleration
        K;                                  % deflection gradient (V/m)
        phik;                               % max deflection phase
        fields;                             % object storing Ex, Ey, and Hz fields
        % simulation parameters
        verbose = true;                     % whether to display progress messages.  Keep on for debugging.  Turn off for doing scans.
        movie_rate = 10;                    % how fast movie plays, make larger for faster
    end
    
    methods
        function obj = simulation
            % Constructs a simulation object.  Sets some of the obvious
            % parameters
            if (obj.verbose) display('initializing simulation ...'); end
            obj.w0 = 2*pi*obj.c0/obj.lambda;
            obj.L = obj.lambda*obj.beta*obj.num_periods;
            obj.gamma = 1/sqrt(1-obj.beta^2);
            addpath(genpath('.'));
        end
        function obj = compute_fields(obj)
            % Computes the fields with an FDFD simulation
            if (obj.verbose) display('computing the field with FDFD ...'); end        
            compute_fields_OO(obj);
        end
        function optimize_radius(obj)
            % scans radius and updates simulation with resulting radius and
            NR = 20;
            Rs = linspace(10e-9, obj.lambda*obj.beta, NR);
            Gs = zeros(NR,1);
            upd = textprogressbar(NR);
            for i = (1:NR)
                upd(i);
                obj.r_pillars = Rs(i);
                obj.compute_fields;                
                Gs(i) = obj.G;
            end
            [~,I] = max(Gs);
            obj.r_pillars = Rs(I);
            obj.compute_fields;
        end
        function [out, trajectory] = propagate_particle(obj, in, phi)
            % Propagates a particle with phase space vector 'in' =
            % [x0,y0,z0, px0, py0, pz0] to give 'out' vector after full
            % propagation with starting electric fields of phase 'phi'            
            if (obj.verbose) display('propagating the field with FDFD ...');  end            
            if (isempty(obj.fields))
                error('need to solve_fields before calling this method')
            end            
            [out, trajectory] = propagate_particle_OO(obj, in, phi);            
        end
        function [outs] = propagate_particle_vectorized(obj, ins)
            % Propagates a particle with phase space vector 'in' =
            % [x0,y0,z0, px0, py0, pz0] to give 'out' vector after full
            % propagation with starting electric fields of phase 'phi'            
            if (obj.verbose) display('propagating the field with FDFD ...');  end            
            if (isempty(obj.fields))
                error('need to solve_fields before calling this method')
            end            
            [outs] = propagate_particle_OO_vectorized(obj, ins);            
        end        
        function obj = make_video(obj,val)
            % Makes a video of the fields specified in argument (Ex, Ey,
            % Ez) as a string
            if (isempty(obj.fields))
                error('need to solve_fields before calling this method')
            end
            if strcmp(val,'Ez')
                field_array = obj.fields.Ez;
            elseif strcmp(val,'Ey')
                field_array = obj.fields.Ey;                
            elseif strcmp(val,'Hx')
                field_array = obj.fields.Hx;                                
            else
                error('need to specify one of Ez, Ey or Hx as an argument to this method')
            end
            scale = max(abs(field_array(:)));
            scale_factor = 0.2;    
            figure(); clf; colormap(redblue); close all; colormap(redblue)
            for i = (1:round(1000*2*pi/obj.movie_rate))
                E = real(field_array*exp(-1i*i*obj.movie_rate/100));
                imagesc(flipud(transpose([E;E;E;E;E;E;E;E])),[-scale*scale_factor,scale*scale_factor]);   
                pause(0.01); clf;
            end
        end        
    end
end