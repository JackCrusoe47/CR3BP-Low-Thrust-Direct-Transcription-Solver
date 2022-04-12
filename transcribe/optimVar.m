function z = optimVar(transcribe,problem)

% - extract parameters
nState = problem.nState;
nControl = problem.nControl;
nSegment = transcribe.param.nSegment;
nOdd = transcribe.param.nOdd;
nNodes = length(transcribe.traj.nodes.t);
timeFixed = problem.flag.timeFixed;

% optim var structure : z = [ state nodes; segment control; tof ]

% segment node states
zState = reshape(transcribe.traj.nodes.x,nState*nNodes,1);
% segment control
zControl = reshape(  transcribe.traj.nodes.u(:, 1:nOdd-1:nSegment*(nOdd-1)), nControl*nSegment, 1 );
% time of flight
if ~timeFixed
    zToF = transcribe.traj.nodes.t(end);
end

% full optim var
z = [ zState; zControl];
if ~timeFixed
    z = [ z; zToF ];
end

end