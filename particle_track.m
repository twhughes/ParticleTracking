% code for testing particle tracker
n_periods = 40;         % number of DLA periods to propagate through
n_per_lam = 300;        % number of grid points per wavelength in FDFD
lambda = 2e-6;          % free space wavelength (m)
radius = 400e-9;        % pillar radius (m)
gap = 400e-9;           % gap size (m)
beta = 0.5;             % electron speed / speed of light
c0 = 3e8;               % speed of light
me = 9e-31;             % mass of electron (kg)
q = 1.6e-19;            % electric charge (C)
E0 = 1e7;               % incident electric field strength (V/m)
phi = pi/3;             % light phase
dl = lambda/n_per_lam;  % grid resolution (m)

% first compute the fields from FDFD (comment this out after the first time
% or else it will slow things down a lot if you re-solve
%[G,phi,fields_out] = compute_fields(radius, gap, beta, n_per_lam);
[Nx,Ny] = size(fields.Ex);

% input phase space vector [x0, y0, px0, py0] (m, m, kg m/s kg m/s)
% NOTE, x0 must be > 0 (dl works)
in = [dl dl*Ny/2-40*dl me*beta*c0 me*beta*c0/100];

%propagate for a bunch of phases
Nphi = 100;
clf; upd = textprogressbar(Nphi);
phis = linspace(0,2*pi,Nphi);
for i = (1:Nphi)
    upd(i);
    phi = phis(i);
    [out,traj] =  propagate_particle(in, fields, n_periods, 0.5, n_per_lam, 1e5, phi);
    hold all;
    plot(traj(1,:), traj(2,:));
end
xlim([0 lambda*beta*n_periods]);
ylim([0 Ny*dl]);
xlabel('x position (m)')
xlabel('y position (m)')
title('trajectories at differing light phase');
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')