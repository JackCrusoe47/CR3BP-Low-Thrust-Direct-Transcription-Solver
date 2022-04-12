function dx = dynamics_CR3BP_Thrust(t,x,u,paramProb)

mu = paramProb.mu;
tstr = paramProb.tstr;
lstr = paramProb.lstr;
Isp = paramProb.Isp;

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

% - dynamics
dx = [
    vx
    vy
    vz
    rx + 2.*vy + ((uNunit.*uThrust.*(ry.*vz - rz.*vy))./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (uThrust.*uWunit.*((vy.*(rx.*vy - ry.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vz.*(rx.*vz - rz.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc) + (uTunit.*uThrust.*vx)./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2)))./astr - (mu.*(mu + rx - 1))./((mu + rx - 1).^2 + ry.^2 + rz.^2).^(3./2) + ((mu + rx).*(mu - 1))./((mu + rx).^2 + ry.^2 + rz.^2).^(3./2)
    ry - 2.*vx + ((uThrust.*uWunit.*((vx.*(rx.*vy - ry.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) - (vz.*(ry.*vz - rz.*vy))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc) - (uNunit.*uThrust.*(rx.*vz - rz.*vx))./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (uTunit.*uThrust.*vy)./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2)))./astr - (mu.*ry)./((mu + rx - 1).^2 + ry.^2 + rz.^2).^(3./2) + (ry.*(mu - 1))./((mu + rx).^2 + ry.^2 + rz.^2).^(3./2)
    ((uThrust.*uWunit.*((vx.*(rx.*vz - rz.*vx))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (vy.*(ry.*vz - rz.*vy))./(((vx).^2 + (vy).^2 + (vz).^2).^(1./2).*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2))))./(1000.*Msc) + (uNunit.*uThrust.*(rx.*vy - ry.*vx))./(1000.*Msc.*((rx.*vy - ry.*vx).^2 + (rx.*vz - rz.*vx).^2 + (ry.*vz - rz.*vy).^2).^(1./2)) + (uTunit.*uThrust.*vz)./(1000.*Msc.*((vx).^2 + (vy).^2 + (vz).^2).^(1./2)))./astr - (mu.*rz)./((mu + rx - 1).^2 + ry.^2 + rz.^2).^(3./2) + (rz.*(mu - 1))./((mu + rx).^2 + ry.^2 + rz.^2).^(3./2)
    -(tstr.*uThrust)./(Isp(uThrust).*g0)
    ];

end