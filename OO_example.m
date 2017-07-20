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
in = [0 0 1.1*p0 0];
% scan through different phis and draw the trajectoroes
NP = 10;
phis = 2*pi*rand(NP,1);
upd = textprogressbar(NP);
s.verbose = false;                  % turn of verbosity (dont print)
clf; hold all;
out_distributions = zeros(4,NP);
for i = (1:NP)
    upd(i);
    [out, traj] = s.propagate_particle(in, phis(i));    
    plot(1e6*traj(1,:),1e6*traj(2,:))
    out_distributions(:,i) = out;
end

xlabel('x position (um)');
ylabel('y position (um)');
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
plot([0,s.Nx*s.dl*s.num_periods*1e6],[-s.gap/2*1e6 -s.gap/2*1e6],'k');
plot([0,s.Nx*s.dl*s.num_periods*1e6],[ s.gap/2*1e6  s.gap/2*1e6],'k');

%%
figure(2); clf;
subplot(3,1,1)
scatter(phis/2/pi,out_distributions(2,:)*1e6,'filled'); 
hold all;
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
ylabel('final y position (um)');
title('transverse position')

subplot(3,1,2)
scatter(phis/2/pi,out_distributions(3,:)/p0,'filled')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
ylabel('x momentum (me \beta c_0)')
title('longitudinal momentum')
subplot(3,1,3)
scatter(phis/2/pi,out_distributions(4,:)/p0,'filled')
xlabel('light phase (2\pi)')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
ylabel('y momentum (me \beta c_0)')
title('transverse momentum')
