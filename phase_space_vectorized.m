% initialize simulation if hasn't already been initialized (clear s)
if (~exist('s'))
    s = simulation;
end

% compute the fields if haven't already been calculated, these are now stored in s along with other goodies
if isempty(s.fields)
    s.compute_fields;
end

p0 = s.gamma*s.me*s.beta*s.c0;
% define input phase space vector
% scan through different phis and draw the trajectoroes
Ne = 5000;                           % number of electrons
s.verbose = false;                  % turn of verbosity (dont print)
clf; hold all;
out_distributions = zeros(4,Ne);
in_distributions = zeros(5,Ne);
subplot(2,2,1); hold all;

E_spread_eV = 10;                           % longitudinal energy spread in eV
p_spread = s.beta*E_spread_eV*s.q/s.c0;     % conversion into momentum spread

ins = zeros(5,Ne);

phi_vec = rand(1,Ne)*2*pi;                              % random phase
pz0_vec = p0*ones(1,Ne) + p_spread*randn(1,Ne);
y0_vec = rand(1,Ne)*s.gap-ones(1,Ne)*s.gap/2;
ins(2,:) = y0_vec;
ins(3,:) = pz0_vec;
ins(5,:) = phi_vec;

outs = s.propagate_particle_vectorized(ins);    

%%
figure(1); clf;
colormap(hsv);
colorbar();


subplot(3,2,1);
scatter(ins(3,:)./p0,outs(4,:)./p0,6,phi_vec,'filled');
xlabel('input momentum (p0)')
ylabel('final transverse momentum (p0)')

subplot(3,2,2);
scatter(ins(3,:)./p0,(outs(3,:)-ins(3,:))./p0,6,phi_vec,'filled');
xlabel('input momentum (p0)')
ylabel('change in longitudinal momentum (p0)')

subplot(3,2,3);
scatter(ins(3,:)./p0,outs(2,:)*1e6,6,phi_vec,'filled');
xlabel('input momentum (p0)')
ylabel('final transverse position (um)')

subplot(3,2,4);
scatter(ins(2,:)*1e6,outs(4,:)./p0,6,phi_vec,'filled');
xlabel('input transverse position (um)')
ylabel('final transverse momentum (p0)')

subplot(3,2,5);
scatter(ins(2,:)*1e6,(outs(3,:)-ins(3,:))./p0,6,phi_vec,'filled');
xlabel('input transverse position (um)')
ylabel('change in longitudinal momentum (p0)')

subplot(3,2,6);
scatter(ins(2,:)*1e6,outs(2,:)*1e6,6,phi_vec,'filled');
xlabel('input transverse position (um)')
ylabel('final transverse position (um)')

