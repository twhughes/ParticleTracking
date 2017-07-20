% initialize simulation if hasn't already been initialized
if (~exist('s'))
    s = simulation;
end

% compute the fields if haven't already been calculated, these are now stored in s along with other goodies
if isempty(s.fields)
    s.compute_fields;
end

p0 = s.beta*s.me*s.c0;
% define input phase space vector
% scan through different phis and draw the trajectoroes
Ne = 1000;                           % number of electrons
upd = textprogressbar(Ne);
s.verbose = false;                  % turn of verbosity (dont print)
clf; hold all;
out_distributions = zeros(4,Ne);
in_distributions = zeros(5,Ne);
subplot(2,2,1); hold all;

for i = (1:Ne)
    upd(i);
    phi = rand()*2*pi;                              % random phase
    x0 = 0;
    y0 = 0;
    px0 = p0 + p0/50*randn();
    py0 = 0;
    in = [x0, y0, px0, py0];
    in_distributions(:,i) = ([in phi]);
    [out, traj] = s.propagate_particle(in, phi);    
    plot(1e6*traj(1,:),1e6*traj(2,:))
    out_distributions(:,i) = out;
end

%%
subplot(2,2,1);
xlabel('x position (um)');
ylabel('y position (um)');
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
plot([0,s.Nx*s.dl*s.num_periods*1e6],[-s.gap/2*1e6 -s.gap/2*1e6],'k');
plot([0,s.Nx*s.dl*s.num_periods*1e6],[ s.gap/2*1e6  s.gap/2*1e6],'k');

subplot(2,2,2)
scatter(in_distributions(3,:)/p0,out_distributions(2,:)*1e6,6,'k','filled'); 
hold all;
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
ylabel('final y position (um)');
xlabel('input momentum (p0)');
title('transverse position')

subplot(2,2,3)
scatter(in_distributions(3,:)/p0,(out_distributions(3,:)-in_distributions(3,:))/p0,6,'k','filled'); 
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
xlabel('input momentum (p0)');
ylabel('momentum change (p0)');
title('change in longitudinal momentum')

subplot(2,2,4)
scatter(in_distributions(3,:)/p0,out_distributions(4,:)/p0,6,'k','filled'); 
xlabel('light phase (2\pi)')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
xlabel('input momentum (p0)');
ylabel('final transverse momentum (p0)');
title('transverse momentum')
%}