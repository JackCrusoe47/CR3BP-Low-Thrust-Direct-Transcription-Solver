function solution = interpTraj(solution,transcribe,problem,solve_options)

% - extract parameters
nOdd = transcribe.param.nOdd;
nOrder = transcribe.param.nOrder;
nSegment = transcribe.param.nSegment;
dtInterp = solve_options.dtInterp;

% - extract state, control and tof
xdata = solution.traj.nodes.x;
udata = solution.traj.nodes.u;
tdata = solution.traj.nodes.t;
tof = tdata(end);

% - segment length
Deltat = tof/transcribe.param.nSegment;

% - number of interpolation points
nInterp = ceil(Deltat/dtInterp);

% - trascription matrices and vectors
Ainv = transcribe.param.Ainv;

% - function handles
dyn = problem.func.dynamics;

% - interpolation points (-1,1)
tauInterp = linspace(-1,1,nInterp);

% - create new B matrix
B = zeros(nOrder+1,nInterp-1);
B(1,:) = 1;
for i = 1:nOrder
    B(1+i,:) = tauInterp(2:nInterp).^(i);
end

xout = xdata(:,1);
uout = udata(:,1);
tout = 0;

idx0_xNodes = 1;
for i = 1:nSegment
    idxf_xNodes = idx0_xNodes + nOdd - 1;
    xNodes = xdata(:,idx0_xNodes:idxf_xNodes);
    
    X = [ xNodes, Deltat/2*dyn(ones(1,nOdd), xNodes, udata(:,i)*ones(1,nOdd) ) ];
    
    C =  X*Ainv;
    
    xout = [xout, C*B];
    uout = [uout, udata(:,i).*ones(1,nInterp-1)];
    tout = [tout, tout(end)+(tauInterp(2:nInterp) + 1)*Deltat/2];
    
    idx0_xNodes = idxf_xNodes;
end

solution.traj.data.x = xout;
solution.traj.data.u = uout;
solution.traj.data.t = tout;
solution.traj.interp.x = @(z) interp1(tout',xout',z','pchip')';
solution.traj.interp.u = @(z) interp1(tout',uout',z','pchip')';

end