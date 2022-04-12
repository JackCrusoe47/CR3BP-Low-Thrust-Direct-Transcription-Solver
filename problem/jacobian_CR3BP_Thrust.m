function dxGrad = jacobian_CR3BP_Thrust(t,x,u,paramProb)

mu = paramProb.mu;
tstr = paramProb.tstr;
lstr = paramProb.lstr;
Isp = paramProb.Isp;
dIsp = paramProb.dIsp;

% acceleration due to gravity
g0 = 9.80665;

% characteristic velocity
vstr = lstr/tstr;

% characteristic acceleration
astr = vstr/tstr;

% number of points
N = length(t);

% - states
rx = x(1,:);
ry = x(2,:);
rz = x(3,:);
vx = x(4,:);
vy = x(5,:);
vz = x(6,:);
Msc = x(7,:);

% - control
uThrust = u(1,:);    % thrust
uTunit = u(2,:);    %
uWunit = u(3,:);
uNunit = u(4,:);

% uThrust = abs(uThrust);

% % - normalizing (for numerical accuracy)
% unorm = vecnorm(u(2:4,:));
% uTunit = uTunit./unorm;
% uWunit = uWunit./unorm;
% uNunit = uNunit./unorm;

% - analytic gradients

% - jacobian dimension
% ( [7 states] x [7 states + 4 control] x [N] )

% - intializing matrix
dxGrad = zeros(7,11,N);

% - derivative wrt state rx
dxGrad(4,1,:) = (mu - 1)./((mu + rx).^2 + ry.^2 + rz.^2).^(3./2) - ((uThrust.*uWunit.*(vy.^2./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + vz.^2./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (vy.*(2.*vy.*(rx.*vy - ry.*vx) + 2.*vz.*(rx.*vz - rz.*vx)).*(rx.*vy - ry.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (vz.*(2.*vy.*(rx.*vy - ry.*vx) + 2.*vz.*(rx.*vz - rz.*vx)).*(rx.*vz - rz.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2))))./(1000.*Msc) + (uNunit.*uThrust.*(2.*vy.*(rx.*vy - ry.*vx) + 2.*vz.*(rx.*vz - rz.*vx)).*(ry.*vz - rz.*vy))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)))./astr - mu./((mu + rx - 1).^2 + ry.^2 + rz.^2).^(3./2) + (3.*mu.*(mu + rx - 1).*(2.*mu + 2.*rx - 2))./(2.*((mu + rx - 1).^2 + ry.^2 + rz.^2).^(5./2)) - (3.*(mu + rx).*(mu - 1).*(2.*mu + 2.*rx))./(2.*((mu + rx).^2 + ry.^2 + rz.^2).^(5./2)) + 1;
dxGrad(5,1,:) = ((uThrust.*uWunit.*((vx.*vy)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (vx.*(2.*vy.*(rx.*vy - ry.*vx) + 2.*vz.*(rx.*vz - rz.*vx)).*(rx.*vy - ry.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (vz.*(2.*vy.*(rx.*vy - ry.*vx) + 2.*vz.*(rx.*vz - rz.*vx)).*(ry.*vz - rz.*vy))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2))))./(1000.*Msc) - (uNunit.*uThrust.*vz)./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (uNunit.*uThrust.*(2.*vy.*(rx.*vy - ry.*vx) + 2.*vz.*(rx.*vz - rz.*vx)).*(rx.*vz - rz.*vx))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)))./astr + (3.*mu.*ry.*(2.*mu + 2.*rx - 2))./(2.*((mu + rx - 1).^2 + ry.^2 + rz.^2).^(5./2)) - (3.*ry.*(mu - 1).*(2.*mu + 2.*rx))./(2.*((mu + rx).^2 + ry.^2 + rz.^2).^(5./2));
dxGrad(6,1,:) = (3.*mu.*rz.*(2.*mu + 2.*rx - 2))./(2.*((mu + rx - 1).^2 + ry.^2 + rz.^2).^(5./2)) - ((uThrust.*uWunit.*((vx.*(2.*vy.*(rx.*vy - ry.*vx) + 2.*vz.*(rx.*vz - rz.*vx)).*(rx.*vz - rz.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (vx.*vz)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vy.*(2.*vy.*(rx.*vy - ry.*vx) + 2.*vz.*(rx.*vz - rz.*vx)).*(ry.*vz - rz.*vy))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2))))./(1000.*Msc) - (uNunit.*uThrust.*vy)./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (uNunit.*uThrust.*(2.*vy.*(rx.*vy - ry.*vx) + 2.*vz.*(rx.*vz - rz.*vx)).*(rx.*vy - ry.*vx))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)))./astr - (3.*rz.*(mu - 1).*(2.*mu + 2.*rx))./(2.*((mu + rx).^2 + ry.^2 + rz.^2).^(5./2));

% - derivative wrt state ry
dxGrad(4,2,:) = ((uNunit.*uThrust.*vz)./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (uThrust.*uWunit.*((vy.*(2.*vx.*(rx.*vy - ry.*vx) - 2.*vz.*(ry.*vz - rz.*vy)).*(rx.*vy - ry.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (vx.*vy)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vz.*(2.*vx.*(rx.*vy - ry.*vx) - 2.*vz.*(ry.*vz - rz.*vy)).*(rx.*vz - rz.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2))))./(1000.*Msc) + (uNunit.*uThrust.*(2.*vx.*(rx.*vy - ry.*vx) - 2.*vz.*(ry.*vz - rz.*vy)).*(ry.*vz - rz.*vy))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)))./astr + (3.*mu.*ry.*(mu + rx - 1))./((mu + rx - 1).^2 + ry.^2 + rz.^2).^(5./2) - (3.*ry.*(mu + rx).*(mu - 1))./((mu + rx).^2 + ry.^2 + rz.^2).^(5./2);
dxGrad(5,2,:) = (mu - 1)./((mu + rx).^2 + ry.^2 + rz.^2).^(3./2) - ((uThrust.*uWunit.*(vx.^2./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + vz.^2./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (vx.*(2.*vx.*(rx.*vy - ry.*vx) - 2.*vz.*(ry.*vz - rz.*vy)).*(rx.*vy - ry.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (vz.*(2.*vx.*(rx.*vy - ry.*vx) - 2.*vz.*(ry.*vz - rz.*vy)).*(ry.*vz - rz.*vy))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2))))./(1000.*Msc) + (uNunit.*uThrust.*(2.*vx.*(rx.*vy - ry.*vx) - 2.*vz.*(ry.*vz - rz.*vy)).*(rx.*vz - rz.*vx))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)))./astr - mu./((mu + rx - 1).^2 + ry.^2 + rz.^2).^(3./2) + (3.*mu.*ry.^2)./((mu + rx - 1).^2 + ry.^2 + rz.^2).^(5./2) - (3.*ry.^2.*(mu - 1))./((mu + rx).^2 + ry.^2 + rz.^2).^(5./2) + 1;
dxGrad(6,2,:) = ((uThrust.*uWunit.*((vy.*vz)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vx.*(2.*vx.*(rx.*vy - ry.*vx) - 2.*vz.*(ry.*vz - rz.*vy)).*(rx.*vz - rz.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (vy.*(2.*vx.*(rx.*vy - ry.*vx) - 2.*vz.*(ry.*vz - rz.*vy)).*(ry.*vz - rz.*vy))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2))))./(1000.*Msc) - (uNunit.*uThrust.*vx)./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (uNunit.*uThrust.*(2.*vx.*(rx.*vy - ry.*vx) - 2.*vz.*(ry.*vz - rz.*vy)).*(rx.*vy - ry.*vx))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)))./astr - (3.*ry.*rz.*(mu - 1))./((mu + rx).^2 + ry.^2 + rz.^2).^(5./2) + (3.*mu.*ry.*rz)./((mu + rx - 1).^2 + ry.^2 + rz.^2).^(5./2);

% - derivative wrt state rz
dxGrad(4,3,:) = (3.*mu.*rz.*(mu + rx - 1))./((mu + rx - 1).^2 + ry.^2 + rz.^2).^(5./2) - ((uThrust.*uWunit.*((vy.*(2.*vx.*(rx.*vz - rz.*vx) + 2.*vy.*(ry.*vz - rz.*vy)).*(rx.*vy - ry.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (vx.*vz)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vz.*(2.*vx.*(rx.*vz - rz.*vx) + 2.*vy.*(ry.*vz - rz.*vy)).*(rx.*vz - rz.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2))))./(1000.*Msc) + (uNunit.*uThrust.*vy)./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (uNunit.*uThrust.*(2.*vx.*(rx.*vz - rz.*vx) + 2.*vy.*(ry.*vz - rz.*vy)).*(ry.*vz - rz.*vy))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)))./astr - (3.*rz.*(mu + rx).*(mu - 1))./((mu + rx).^2 + ry.^2 + rz.^2).^(5./2);
dxGrad(5,3,:) = ((uThrust.*uWunit.*((vy.*vz)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vx.*(2.*vx.*(rx.*vz - rz.*vx) + 2.*vy.*(ry.*vz - rz.*vy)).*(rx.*vy - ry.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (vz.*(2.*vx.*(rx.*vz - rz.*vx) + 2.*vy.*(ry.*vz - rz.*vy)).*(ry.*vz - rz.*vy))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2))))./(1000.*Msc) + (uNunit.*uThrust.*vx)./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (uNunit.*uThrust.*(2.*vx.*(rx.*vz - rz.*vx) + 2.*vy.*(ry.*vz - rz.*vy)).*(rx.*vz - rz.*vx))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)))./astr - (3.*ry.*rz.*(mu - 1))./((mu + rx).^2 + ry.^2 + rz.^2).^(5./2) + (3.*mu.*ry.*rz)./((mu + rx - 1).^2 + ry.^2 + rz.^2).^(5./2);
dxGrad(6,3,:) = (mu - 1)./((mu + rx).^2 + ry.^2 + rz.^2).^(3./2) - ((uThrust.*uWunit.*(vx.^2./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + vy.^2./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (vx.*(2.*vx.*(rx.*vz - rz.*vx) + 2.*vy.*(ry.*vz - rz.*vy)).*(rx.*vz - rz.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (vy.*(2.*vx.*(rx.*vz - rz.*vx) + 2.*vy.*(ry.*vz - rz.*vy)).*(ry.*vz - rz.*vy))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2))))./(1000.*Msc) - (uNunit.*uThrust.*(2.*vx.*(rx.*vz - rz.*vx) + 2.*vy.*(ry.*vz - rz.*vy)).*(rx.*vy - ry.*vx))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)))./astr - mu./((mu + rx - 1).^2 + ry.^2 + rz.^2).^(3./2) + (3.*mu.*rz.^2)./((mu + rx - 1).^2 + ry.^2 + rz.^2).^(5./2) - (3.*rz.^2.*(mu - 1))./((mu + rx).^2 + ry.^2 + rz.^2).^(5./2);

% - derivative wrt state vx
dxGrad(1,4,:) = ones([1,1,N]);
dxGrad(4,4,:) = ((uTunit.*uThrust)./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2)) + (uThrust.*uWunit.*((ry.*vy)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (rz.*vz)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (vy.*(2.*ry.*(rx.*vy - ry.*vx) + 2.*rz.*(rx.*vz - rz.*vx)).*(rx.*vy - ry.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (vz.*(2.*ry.*(rx.*vy - ry.*vx) + 2.*rz.*(rx.*vz - rz.*vx)).*(rx.*vz - rz.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (vy.*(vx).*(rx.*vy - ry.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vz.*(vx).*(rx.*vz - rz.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc) + (uNunit.*uThrust.*(2.*ry.*(rx.*vy - ry.*vx) + 2.*rz.*(rx.*vz - rz.*vx)).*(ry.*vz - rz.*vy))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (uTunit.*uThrust.*vx.*(vx))./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(3./2)))./astr;
dxGrad(5,4,:) = ((uThrust.*uWunit.*((rx.*vy - ry.*vx)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (ry.*vx)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vx.*(2.*ry.*(rx.*vy - ry.*vx) + 2.*rz.*(rx.*vz - rz.*vx)).*(rx.*vy - ry.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (vz.*(2.*ry.*(rx.*vy - ry.*vx) + 2.*rz.*(rx.*vz - rz.*vx)).*(ry.*vz - rz.*vy))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (vx.*(vx).*(rx.*vy - ry.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vz.*(vx).*(ry.*vz - rz.*vy))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc) + (rz.*uNunit.*uThrust)./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (uNunit.*uThrust.*(2.*ry.*(rx.*vy - ry.*vx) + 2.*rz.*(rx.*vz - rz.*vx)).*(rx.*vz - rz.*vx))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (uTunit.*uThrust.*vy.*(vx))./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(3./2)))./astr - 2;
dxGrad(6,4,:) = ((uThrust.*uWunit.*((rx.*vz - rz.*vx)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (rz.*vx)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vx.*(2.*ry.*(rx.*vy - ry.*vx) + 2.*rz.*(rx.*vz - rz.*vx)).*(rx.*vz - rz.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (vy.*(2.*ry.*(rx.*vy - ry.*vx) + 2.*rz.*(rx.*vz - rz.*vx)).*(ry.*vz - rz.*vy))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (vx.*(vx).*(rx.*vz - rz.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (vy.*(vx).*(ry.*vz - rz.*vy))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc) - (ry.*uNunit.*uThrust)./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (uNunit.*uThrust.*(2.*ry.*(rx.*vy - ry.*vx) + 2.*rz.*(rx.*vz - rz.*vx)).*(rx.*vy - ry.*vx))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (uTunit.*uThrust.*vz.*(vx))./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(3./2)))./astr;

% - derivative wrt state vy
dxGrad(2,5,:) = ones([1,1,N]);
dxGrad(4,5,:) = 2 - ((rz.*uNunit.*uThrust)./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (uThrust.*uWunit.*((vy.*(2.*rx.*(rx.*vy - ry.*vx) - 2.*rz.*(ry.*vz - rz.*vy)).*(rx.*vy - ry.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (rx.*vy)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (rx.*vy - ry.*vx)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vz.*(2.*rx.*(rx.*vy - ry.*vx) - 2.*rz.*(ry.*vz - rz.*vy)).*(rx.*vz - rz.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (vy.*(vy).*(rx.*vy - ry.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vz.*(vy).*(rx.*vz - rz.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc) + (uNunit.*uThrust.*(2.*rx.*(rx.*vy - ry.*vx) - 2.*rz.*(ry.*vz - rz.*vy)).*(ry.*vz - rz.*vy))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (uTunit.*uThrust.*vx.*(vy))./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(3./2)))./astr;
dxGrad(5,5,:) = ((uTunit.*uThrust)./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2)) + (uThrust.*uWunit.*((rx.*vx)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (rz.*vz)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (vx.*(2.*rx.*(rx.*vy - ry.*vx) - 2.*rz.*(ry.*vz - rz.*vy)).*(rx.*vy - ry.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (vz.*(2.*rx.*(rx.*vy - ry.*vx) - 2.*rz.*(ry.*vz - rz.*vy)).*(ry.*vz - rz.*vy))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (vx.*(vy).*(rx.*vy - ry.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vz.*(vy).*(ry.*vz - rz.*vy))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc) + (uNunit.*uThrust.*(2.*rx.*(rx.*vy - ry.*vx) - 2.*rz.*(ry.*vz - rz.*vy)).*(rx.*vz - rz.*vx))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (uTunit.*uThrust.*vy.*(vy))./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(3./2)))./astr;
dxGrad(6,5,:) = -((uThrust.*uWunit.*((rz.*vy)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (ry.*vz - rz.*vy)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vx.*(2.*rx.*(rx.*vy - ry.*vx) - 2.*rz.*(ry.*vz - rz.*vy)).*(rx.*vz - rz.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (vy.*(2.*rx.*(rx.*vy - ry.*vx) - 2.*rz.*(ry.*vz - rz.*vy)).*(ry.*vz - rz.*vy))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (vx.*(vy).*(rx.*vz - rz.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vy.*(vy).*(ry.*vz - rz.*vy))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc) - (rx.*uNunit.*uThrust)./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (uNunit.*uThrust.*(2.*rx.*(rx.*vy - ry.*vx) - 2.*rz.*(ry.*vz - rz.*vy)).*(rx.*vy - ry.*vx))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (uTunit.*uThrust.*vz.*(vy))./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(3./2)))./astr;

% - derivative wrt state vz
dxGrad(3,6,:) = ones([1,1,N]);
dxGrad(4,6,:) = ((uThrust.*uWunit.*((vy.*(2.*rx.*(rx.*vz - rz.*vx) + 2.*ry.*(ry.*vz - rz.*vy)).*(rx.*vy - ry.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (rx.*vz)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (rx.*vz - rz.*vx)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vz.*(2.*rx.*(rx.*vz - rz.*vx) + 2.*ry.*(ry.*vz - rz.*vy)).*(rx.*vz - rz.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (vy.*(vz).*(rx.*vy - ry.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vz.*(vz).*(rx.*vz - rz.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc) + (ry.*uNunit.*uThrust)./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (uNunit.*uThrust.*(2.*rx.*(rx.*vz - rz.*vx) + 2.*ry.*(ry.*vz - rz.*vy)).*(ry.*vz - rz.*vy))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (uTunit.*uThrust.*vx.*(vz))./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(3./2)))./astr;
dxGrad(5,6,:) = -((uThrust.*uWunit.*((ry.*vz - rz.*vy)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (ry.*vz)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vx.*(2.*rx.*(rx.*vz - rz.*vx) + 2.*ry.*(ry.*vz - rz.*vy)).*(rx.*vy - ry.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (vz.*(2.*rx.*(rx.*vz - rz.*vx) + 2.*ry.*(ry.*vz - rz.*vy)).*(ry.*vz - rz.*vy))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (vx.*(vz).*(rx.*vy - ry.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (vz.*(vz).*(ry.*vz - rz.*vy))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc) + (rx.*uNunit.*uThrust)./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (uNunit.*uThrust.*(2.*rx.*(rx.*vz - rz.*vx) + 2.*ry.*(ry.*vz - rz.*vy)).*(rx.*vz - rz.*vx))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (uTunit.*uThrust.*vy.*(vz))./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(3./2)))./astr;
dxGrad(6,6,:) = -((uThrust.*uWunit.*((vx.*(2.*rx.*(rx.*vz - rz.*vx) + 2.*ry.*(ry.*vz - rz.*vy)).*(rx.*vz - rz.*vx))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) - (ry.*vy)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (rx.*vx)./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vy.*(2.*rx.*(rx.*vz - rz.*vx) + 2.*ry.*(ry.*vz - rz.*vy)).*(ry.*vz - rz.*vy))./(2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (vx.*(vz).*(rx.*vz - rz.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vy.*(vz).*(ry.*vz - rz.*vy))./(((vx).^2 + (vy).^2 + (vz).^2).^(3./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc) - (uTunit.*uThrust)./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2)) + (uNunit.*uThrust.*(2.*rx.*(rx.*vz - rz.*vx) + 2.*ry.*(ry.*vz - rz.*vy)).*(rx.*vy - ry.*vx))./(2000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(3./2)) + (uTunit.*uThrust.*vz.*(vz))./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(3./2)))./astr;

% - derivative wrt state Msc

dxGrad(4,7,:) = -((uNunit.*uThrust.*(ry.*vz - rz.*vy))./(1000.*Msc.^2.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (uThrust.*uWunit.*((vy.*(rx.*vy - ry.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vz.*(rx.*vz - rz.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc.^2) + (uTunit.*uThrust.*vx)./(1000.*Msc.^2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2)))./astr;
dxGrad(5,7,:) = -((uThrust.*uWunit.*((vx.*(rx.*vy - ry.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (vz.*(ry.*vz - rz.*vy))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc.^2) - (uNunit.*uThrust.*(rx.*vz - rz.*vx))./(1000.*Msc.^2.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (uTunit.*uThrust.*vy)./(1000.*Msc.^2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2)))./astr;
dxGrad(6,7,:) = -((uThrust.*uWunit.*((vx.*(rx.*vz - rz.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vy.*(ry.*vz - rz.*vy))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc.^2) + (uNunit.*uThrust.*(rx.*vy - ry.*vx))./(1000.*Msc.^2.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (uTunit.*uThrust.*vz)./(1000.*Msc.^2.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2)))./astr;

% - derivative wrt control uThrust
dxGrad(4,8,:) = ((uNunit.*(ry.*vz - rz.*vy))./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (uWunit.*((vy.*(rx.*vy - ry.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vz.*(rx.*vz - rz.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc) + (uTunit.*vx)./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2)))./astr;
dxGrad(5,8,:) = ((uWunit.*((vx.*(rx.*vy - ry.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (vz.*(ry.*vz - rz.*vy))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc) - (uNunit.*(rx.*vz - rz.*vx))./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (uTunit.*vy)./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2)))./astr;
dxGrad(6,8,:) = ((uWunit.*((vx.*(rx.*vz - rz.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vy.*(ry.*vz - rz.*vy))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc) + (uNunit.*(rx.*vy - ry.*vx))./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (uTunit.*vz)./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2)))./astr;
dxGrad(7,8,:) = -tstr./(Isp(uThrust).*g0) + (tstr.*uThrust)./(Isp(uThrust).^2.*g0).*dIsp(uThrust);

% - derivative wrt control uTunit
dxGrad(4,9,:) = (uThrust.*vx)./(1000.*Msc.*astr.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2));
dxGrad(5,9,:) = (uThrust.*vy)./(1000.*Msc.*astr.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2));
dxGrad(6,9,:) = (uThrust.*vz)./(1000.*Msc.*astr.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2));

% - derivative wrt control uWunit
dxGrad(4,10,:) = -(uThrust.*((vy.*(rx.*vy - ry.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vz.*(rx.*vz - rz.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc.*astr);
dxGrad(5,10,:) = (uThrust.*((vx.*(rx.*vy - ry.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (vz.*(ry.*vz - rz.*vy))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc.*astr);
dxGrad(6,10,:) = (uThrust.*((vx.*(rx.*vz - rz.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vy.*(ry.*vz - rz.*vy))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc.*astr);

% - derivative wrt control uNunit
dxGrad(4,11,:) = (uThrust.*(ry.*vz - rz.*vy))./(1000.*Msc.*astr.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2));
dxGrad(5,11,:) = -(uThrust.*(rx.*vz - rz.*vx))./(1000.*Msc.*astr.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2));
dxGrad(6,11,:) = (uThrust.*(rx.*vy - ry.*vx))./(1000.*Msc.*astr.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2));

end