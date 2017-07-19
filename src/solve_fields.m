function [] = solve_fields(obj)
    % function to compute the gradient, phase, and EM fields for a given
    % dual pillar simulation.  The pillars are assumed to be made of Si,
    % the wavelength is 2um.
       
    % simulation variables
    obj.lambda = 2e-6;          % wavelength
    obj.npml = 25;              % number of PML grids ( should be > 15 to be safe)
    eps = 3.45^2;           % permittivity of pillars   
   
    % define relevant lengths
    spc_pml_src = obj.lambda;   % space between pml and current source (m)
    spc_src_pil = obj.lambda;   % space between current source and pillar edge (m)
    pil_radius = obj.r_pillars;    % pillar radius (m)
    spc_gap = obj.gap;          % gap distance (m)

    % compute simulation parameters
    dl = obj.lambda/obj.n_per_lam;  % compute grid resolution
    obj.dl = dl;
    RES = [dl dl];          % load into RES variable for FDFD function
    NPML = [0 0 obj.npml obj.npml]; % load PML into NPML varible for FDFD ([-x +x -y +y])
    BC = [-1 -1];           % periodic boundary conditions
        
    % put these lengths in terms of discrete # of grid points
    pts_pml_src = round(spc_pml_src/dl);
    pts_src_tooth = round(spc_src_pil/dl);
    pts_pil = round(pil_radius/dl);
    pts_gap = round(spc_gap/dl);

    % compute number of grid points needed and corresoponding half grids
    Nx = round(obj.lambda*obj.beta/dl);
    Ny = round(2*obj.npml + 2*pts_pml_src + 2*pts_src_tooth + 4*pts_pil + pts_gap);   
    nx = round(Nx/2);
    ny = round(Ny/2);
    obj.Nx = Nx; obj.Ny = Ny; obj.nx = nx; obj.ny = ny; obj.W = dl*Ny;
    % initialize the perimittivity, permeability and Mz current source
    % grids
    ER = ones(Nx,Ny);
    MuR = ones(Nx,Ny);
    b = zeros(Nx,Ny);

    % inject the current sources
    b(:,obj.npml + pts_pml_src) = 1;
    b(:,end-obj.npml-pts_pml_src+1) = 1;  % comment this line out for single side drive

    % solve for the fields in a vacuum for normalization
    [fields, ~] = FDFD(ER,MuR,RES,NPML,BC,obj.lambda,'Hz',b);    % FDFD simulation 
    E_abs = sqrt(abs(fields.Ex).^2 + abs(fields.Ey).^2);    % compute abs(E) everywhere
    E0 = mean(E_abs(:,ny));                                 % treat E0 as the average along gap

    % left and right y positions of the pillars
    y_pil_left = obj.npml + pts_pml_src + pts_src_tooth + pts_pil;
    y_pil_right = Ny - (obj.npml + pts_pml_src + pts_src_tooth + pts_pil) + 1;
    % draw the pillars ono the epsilon grid
    for xi = (1:Nx)
        for yi = (1:Ny)
            if ((xi-nx)^2 + (yi-y_pil_left)^2 < pts_pil^2)
                ER(xi,yi) = eps;
            end
            if ((xi-nx)^2 + (yi-y_pil_right)^2 < pts_pil^2)
                ER(xi,yi) = eps;        
            end                
        end
    end

    % compute the fields (and normalize by E0) for the real permittivity
    % grid (with dual pillars)
    [fields, ~] = FDFD(ER,MuR,RES,NPML,BC,obj.lambda,'Hz',b);
    Ex = fields.Ex/E0;
    Ey = fields.Ey/E0;
    Hz = fields.Hz/E0;
    
    % make field_out variable to return normalized fields
    obj.fields = {};
    obj.fields.Ex = Ex;
    obj.fields.Ey = Ey;
    obj.fields.Hz = Hz;
    
    % get accelerating fields at gap center, define acceleration 'kernel' e
    Ex_mid = Ex(:,ny);
    eta = exp(1i*2*pi*(1:Nx)/Nx)/Nx;

    % compute gradient phasor
    g = sum(Ex_mid.*transpose(eta));
    % maximum acceleration gradient is the absolute value
    obj.G = abs(g);
    % grab the phase at maximum
    obj.phi0 = angle(g);
    
end