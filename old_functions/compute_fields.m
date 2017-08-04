function [G,phi,fields_out] = compute_fields(radius, gap, beta, n_per_lam)
    % function to compute the gradient, phase, and EM fields for a given
    % dual pillar simulation.  The pillars are assumed to be made of Si,
    % the wavelength is 2um.
       
    % simulation variables
    lambda = 2e-6;          % wavelength
    npml = 25;              % number of PML grids ( should be > 15 to be safe)
    eps = 3.45^2;           % permittivity of pillars   
    Pol = 'Hz';             % polarization ('Hz' or 'Ez')
   
    % define relevant lengths
    spc_pml_src = lambda;   % space between pml and current source (m)
    spc_src_pil = lambda;   % space between current source and pillar edge (m)
    pil_radius = radius;    % pillar radius (m)
    spc_gap = gap;          % gap distance (m)

    % compute simulation parameters
    dl = lambda/n_per_lam;  % compute grid resolution
    RES = [dl dl];          % load into RES variable for FDFD function
    NPML = [0 0 npml npml]; % load PML into NPML varible for FDFD ([-z +z -y +y])
    BC = [-1 -1];           % periodic boundary conditions
        
    % put these lengths in terms of discrete # of grid points
    pts_pml_src = round(spc_pml_src/dl);
    pts_src_tooth = round(spc_src_pil/dl);
    pts_pil = round(pil_radius/dl);
    pts_gap = round(spc_gap/dl);

    % compute number of grid points needed and corresoponding half grids
    Nz = round(lambda*beta/dl);
    Ny = round(2*npml + 2*pts_pml_src + 2*pts_src_tooth + 4*pts_pil + pts_gap);   
    nz = round(Nz/2);
    ny = round(Ny/2);

    % initialize the perimittivity, permeability and Mz current source
    % grids
    ER = ones(Nz,Ny);
    MuR = ones(Nz,Ny);
    b = zeros(Nz,Ny);

    % inject the current sources
    b(:,npml + pts_pml_src) = 1;
    b(:,end-npml-pts_pml_src+1) = 1;  % comment this line out for single side drive

    % solve for the fields in a vacuum for normalization
    [fields, ~] = FDFD(ER,MuR,RES,NPML,BC,lambda,Pol,b);    % FDFD simulation 
    E_abs = sqrt(abs(fields.Ex).^2 + abs(fields.Ey).^2);    % compute abs(E) everywhere
    E0 = mean(E_abs(:,ny));                                 % treat E0 as the average along gap

    % left and right y positions of the pillars
    y_pil_left = npml + pts_pml_src + pts_src_tooth + pts_pil;
    y_pil_right = Ny - (npml + pts_pml_src + pts_src_tooth + pts_pil) + 1;
    % draw the pillars ono the epsilon grid
    for zi = (1:Nz)
        for yi = (1:Ny)
            if ((zi-nz)^2 + (yi-y_pil_left)^2 < pts_pil^2)
                ER(zi,yi) = eps;
            end
            if ((zi-nz)^2 + (yi-y_pil_right)^2 < pts_pil^2)
                ER(zi,yi) = eps;        
            end                
        end
    end

    % compute the fields (and normalize by E0) for the real permittivity
    % grid (with dual pillars)
    [fields, ~] = FDFD(ER,MuR,RES,NPML,BC,lambda,Pol,b);
    Ez = fields.Ex/E0;
    Ey = fields.Ey/E0;
    Hx = fields.Hz/E0;
    
    % make field_out variable to return normalized fields
    % NOTE: in my FDFD solver convention, (x,y) are the dimensions in the
    % plane whereas (z) is out of the plane.  We change this so that the z
    % is now the propagation direction.  THIS IS NOT A TYPO!
    fields_out = {};
    fields_out.Ez = Ez;
    fields_out.Ey = Ey;
    fields_out.Hx = Hx;
    
    % get accelerating fields at gap center, define acceleration 'kernel' e
    Ez_mid = Ez(:,ny);
    eta = exp(1i*2*pi*(1:Nz)/Nz)/Nz;

    % compute gradient phasor
    g = sum(Ez_mid.*transpose(eta));
    % maximum acceleration gradient is the absolute value
    G = abs(g);
    % grab the phase at maximum
    phi = angle(g);
    
end