% This is to check what grid resolution we need to model the dual pillars.
% If the gradient changes with grid resolution, this is not good and we
% need to make it smaller until convergence.  Inspecting the plot, it looks
% like this occurs around 100 grid points per wavelength.  Or a dl of
% 2um/200 = 10 nm.

addpath(genpath('./dependencies'));

N = 50;
upd = textprogressbar(N);

ginls = linspace(20,300,N);
R = 500e-9;
gap = 400e-9;

Gs = zeros(N,1);
for i = (1:N) 
    upd(i);
    n_per_lam = ginls(i);
    [G,phi,fields] = compute_fields(R, gap, beta, n_per_lam);
    Gs(i) = G; 
end

%%

% make a video of the fields
figure(1); clf; colormap(redblue);
for i = (1:1000)
    E = (real(fields.Ex*exp(1i*i/10))); imagesc([E;E;E;E;E;E;E;E]);
    pause(0.01); clf;
end

figure(2); clf;
plot(ginls,Gs);
xlabel('grid points per \lambda');
ylabel('gradient (E_0)');
title('checking convergence of gradient');
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
xlim([min(ginls) max(ginls)]);