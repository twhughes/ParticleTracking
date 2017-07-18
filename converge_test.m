% This is to check what grid resolution we need to model the dual pillars.
% If the gradient changes with grid resolution, this is not good and we
% need to make it smaller until convergence.  Inspecting the plot, it looks
% like this occurs around 100 grid points per wavelength.  Or a dl of
% 2um/200 = 10 nm.

% add all dependencies
addpath(genpath('./dependencies'));

% number of grid resolutions to look at
N = 50;
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

% make a video of the fields in time
figure(1); clf; colormap(redblue);
for i = (1:400)
    E = (real(fields.Ex*exp(-1i*i/10))); 
    imagesc(flipud(transpose([E;E;E;E;E;E;E;E])),[-1.5,1.5]);
    title('E_x');
    xlabel('x grid points');
    ylabel('y grid points');
    colorbar();    
    pause(0.01); clf;
end

% plot the convergence
figure(2); clf;
plot(ginls,Gs);
xlabel('grid points per \lambda');
ylabel('gradient (E_0)');
title('checking convergence of gradient');
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
xlim([min(ginls) max(ginls)]);