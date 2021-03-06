% code for testing particle tracker
n_periods = 5;         % number of DLA periods to propagate through
n_per_lam = 300;        % number of grid points per wavelength in FDFD
lambda = 2e-6;          % free space wavelength (m)
radius = 400e-9;        % pillar radius (m)
gap = 400e-9;           % gap size (m)
beta = 0.5;             % electron speed / speed of light
c0 = 3e8;               % speed of light
me = 9e-31;             % mass of electron (kg)
q = 1.6e-19;            % electric charge (C)
E0 = 1e6;               % incident electric field strength (V/m)
phi = pi/3;             % light phase
dl = lambda/n_per_lam;  % grid resolution (m)

% first compute the fields from FDFD (comment this out after the first time
% or else it will slow things down a lot if you re-solve
%[G,phi0,fields_out] = compute_fields(radius, gap, beta, n_per_lam);
[Nx,Ny] = size(fields.Ex);

% lets choose to work in units of kg <=> me, m <=> dl, sec <=> dl/c0
% input phase space vector [x0, y0, px0, py0] (m, m, kg m/s kg m/s)
% NOTE, x0 must be > 0 (dl works)
in = [dl dl*Ny/2+dl*4 me*beta*c0 0*me*beta*c0];

%propagate for a bunch of phases
Nphi = 1000;
clf; upd = textprogressbar(Nphi);
phis = 2*pi*rand(Nphi,1);
final_pxs = zeros(Nphi,1);
final_ys = zeros(Nphi,1);

for i = (1:Nphi)
    upd(i);
    phi = phis(i)-phi0;
    [out,traj] =  propagate_particle(in, fields, n_periods, beta, n_per_lam, E0, phi);
    hold all;
    plot(traj(1,:), traj(2,:));
    final_pxs(i) = out(3);
    final_ys(i) = out(2);
    
end

xlim([0 lambda*beta*n_periods]);
ylim([0 Ny*dl]);
xlabel('x position (m)')
ylabel('y position (m)')
title('trajectories at differing light phase');
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')