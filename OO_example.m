% initialize simulation if hasn't already been initialized
if (~exist('s'))
    s = simulation;
end

s.num_periods = 100;
% compute the fields if haven't already been calculated, these are now stored in s along with other goodies
if isempty(s.fields)
    s.compute_fields;
end

% define input phase space vector
in = [s.dl s.W/2 s.beta*s.me*s.c0 0];
% scan through different phis and draw the trajectoroes
NP = 5;
phis = 2*pi*rand(NP,1);
upd = textprogressbar(NP);
s.verbose = false;                  % turn of verbosity (dont print)
clf; hold all;
for i = (1:NP)
    upd(i);
    [out, traj] = s.propagate_particle(in, phis(i));    
    plot(1e6*traj(1,:),1e6*traj(2,:))
end

xlabel('x position (um)');
ylabel('y position (um)');
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')