%% cleanup and startup

clear; clc; close all;

addpath(genpath('..'));

%% problem definition


% problem parameters
problem.param.mu = 0.012150584269940;                        % Earth-Moon CR3BP characteristic parameter         
problem.param.tstr = 3.756998590391565e+05;                  % Earth-Moon CR3BP characteristic time  
problem.param.lstr = 3.847479919800000e+05;                  % Earth-Moon CR3BP characteristic length 
problem.param.Isp = @(thrust) 2700 - 300/0.3e-3.*thrust;     % Specific impulse function
problem.param.dIsp = @(thrust) - 300/0.3e-3;                 % Specific impulse derivative wrt. thrust
problem.param.maxThrust = 1.2e-3;                            % maximum spacecraft thrust
problem.param.maxMass = 24;                                  % maximum spacecraft mass 

% CR3BP intial state and time of flight
SV0 = [1.021881345449140;0;-0.182000000000000;0;-0.102950816697240;0];
tf = 1.509263667196279;

% full state of CR3BP+low-thrust problem
x0 = [SV0;problem.param.maxMass];

% control (retrograde burn at 1/4th max thrust)
u0 = [problem.param.maxThrust/4;-1;0;0];

% running ode solver to get initial trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-14);
[t,x] = ode45(@(t,x)dynamics_CR3BP_Thrust(t,x,u0,problem.param),[0,2*tf],x0,options);
x = x';

% adding trajectory to struct
problem.traj.data.t = t;
problem.traj.data.x = x;
problem.traj.data.u = u0*ones([1,length(t)]);

% dimensions of problem
problem.nState = 7;
problem.nControl = 4;

% problem functions
problem.func.dynamics = @(t,x,u)dynamics_CR3BP_Thrust(t,x,u,problem.param);
problem.func.jacobian = @(t,x,u)jacobian_CR3BP_Thrust(t,x,u,problem.param);

% time fixed flag
problem.flag.timeFixed = false;

% new target location
dx = 10000/384600;      % +10000km on x-axis
dy = 2000/384600;       % +200km on y-axis
dz = -1000/384600;      % -1000km on z-axis

% problem bounds
problem.bnds.path0.lb = x(:,1);
problem.bnds.path0.ub = x(:,1);
problem.bnds.pathf.lb = [ x(1,end)+dx; x(2,end)+dy; x(3,end)+dz; -inf; -inf; -inf; 0 ];
problem.bnds.pathf.ub = [ x(1,end)+dx; x(2,end)+dy; x(3,end)+dz; +inf; +inf; +inf; 24 ];

%% solver setup

solver_options.maxIter = 100;
solver_options.tolFunc = 1e-10;
solver_options.dtInterp = 0.001;

%% solve feasible

problem1 = problem;

problem1.objType = 'feasible';
solver_options.nSegment = 100;
solver_options.nOrder = 7;
[solution1,problem1] = optimalTranscribtion(problem1,solver_options);

%% solve optimal

problem2 = problem;
problem2.traj.data = solution1.traj.data;

problem2.objType = 'thrustOptimal';
solver_options.nSegment = 120;
solver_options.nOrder = 11;
solver_options.maxIter = 500;
[solution2,problem2] = optimalTranscribtion(problem2,solver_options);


%% plot solution

% -- update with new feasible trajectory
figure(); hold on; grid on;
% axis equal;
title('Trajectory');
xlabel('x-axis [nd]'); ylabel('y-axis [nd]'); zlabel('z-axis [nd]');
plot3(1-problem.param.mu,0,0,'ok','MarkerFaceColor',[0.7;0.7;0.7],'MarkerSize',10,'HandleVisibility','off');
text(1-problem.param.mu,0,0,'   Moon','HandleVisibility','off');
plot3(problem.traj.data.x(1,1),problem.traj.data.x(2,1),problem.traj.data.x(3,1),'.k','HandleVisibility','off');
text(problem.traj.data.x(1,1),problem.traj.data.x(2,1),problem.traj.data.x(3,1),'  t0','HandleVisibility','off');
plot3(problem.traj.data.x(1,:),problem.traj.data.x(2,:),problem.traj.data.x(3,:),'--k');
plot3(problem.traj.data.x(1,end),problem.traj.data.x(2,end),problem.traj.data.x(3,end),'*k','HandleVisibility','off');
text(problem.traj.data.x(1,end),problem.traj.data.x(2,end),problem.traj.data.x(3,end),'  tf','HandleVisibility','off');
plot3(solution1.traj.data.x(1,end),solution1.traj.data.x(2,end),solution1.traj.data.x(3,end),'*g','HandleVisibility','off');
text(solution1.traj.data.x(1,end),solution1.traj.data.x(2,end),solution1.traj.data.x(3,end),' target','Color','g','HandleVisibility','off');
plot3(solution1.traj.data.x(1,:),solution1.traj.data.x(2,:),solution1.traj.data.x(3,:),'b');
plot3(solution2.traj.data.x(1,:),solution2.traj.data.x(2,:),solution2.traj.data.x(3,:),'r');
legend('guess','feasible','optimized');
view(30,30);

figure();
subplot(2,1,1);
hold on; grid on;title('Mass');
xlabel('time [nd]');
ylabel('mass [kg]');
plot(problem.traj.data.t,problem.traj.data.x(7,:),'k');
plot(solution1.traj.data.t,solution1.traj.data.x(7,:),'b');
plot(solution2.traj.data.t,solution2.traj.data.x(7,:),'r');
legend('guess','feasible','optimized');

subplot(2,1,2);
xlabel('time [nd]');
ylabel('thrust [N]');
hold on; grid on;title('Thrust');
plot(problem.traj.data.t,problem.traj.data.u(1,:),'k');
plot(solution1.traj.data.t,solution1.traj.data.u(1,:),'b');
plot(solution2.traj.data.t,solution2.traj.data.u(1,:),'r');
legend('guess','feasible','optimized');
