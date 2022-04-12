function [z_lb,z_ub] = bndsVar(transcribe,problem)

% - extract parameters
nState = problem.nState;
nControl = problem.nControl;
nSegment = transcribe.param.nSegment;
% nOdd = transcribe.param.nOdd;
nNodes = length(transcribe.traj.nodes.t);
timeFixed = problem.flag.timeFixed;

% optim var structure : z = [ state nodes; segment control; tof ]

x.lb = problem.bnds.state.lb;
x.ub = problem.bnds.state.ub;
u.lb = problem.bnds.control.lb;
u.ub = problem.bnds.control.ub;
if ~timeFixed
    tf.lb = problem.bnds.time.lb;
    tf.ub = problem.bnds.time.ub;
end

z_lb = [ reshape(x.lb.*ones(1,nNodes),nState*nNodes,1);
    reshape(u.lb.*ones(1,nSegment),nControl*nSegment,1) ];

z_ub = [ reshape(x.ub.*ones(1,nNodes),nState*nNodes,1);
    reshape(u.ub.*ones(1,nSegment),nControl*nSegment,1) ];

if ~timeFixed
    z_lb = [ z_lb; tf.lb];
    z_ub = [ z_ub; tf.ub];
end

end