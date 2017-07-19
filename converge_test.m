% This is to check what grid resolution we need to model the dual pillars.
% If the gradient changes with grid resolution, this is not good and we
% need to make it smaller until convergence.  Inspecting the plot, it looks
% like this occurs around 100 grid points per wavelength.  Or a dl of
% 2um/200 = 10 nm.

% add all dependencies
addpath(genpath('./dependencies'));

% number of grid resolutions to look at
N = 10;
% fancy text progress bar (ignore)
upd = textprogressbar(N);

% grids per wavelength to scan
ginls = linspace(20,300,N);

R = 500e-9;     % pillar radius (m)
gap = 400e-9;   % gap spaing (m)
beta = 0.5;     % electron speed / speed of light

Gs = zeros(N,1);  % store gradients
for i = (1:N) 
    upd(i);                     % update progress bar
    n_per_lam = ginls(i);       % get resolution
    [G,phi,fields] = compute_fields(R, gap, beta, n_per_lam);  % FDFD simulation
    Gs(i) = G;                  % save gradient (in units of E0)
end

%%

% plot the convergence
figure(1); clf; close all;
plot(ginls,Gs);
xlabel('grid points per \lambda');
ylabel('gradient (E_0)');
title('checking convergence of gradient');
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
xlim([min(ginls) max(ginls)]);

make_video(fields.Ex);
make_video(fields.Ey);
make_video(fields.Hz);