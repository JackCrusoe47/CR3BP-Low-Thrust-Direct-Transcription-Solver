function [F_lb,F_ub] = bndsConst(transcribe,problem)

% - extract parameters
nState = problem.nState;
% nControl = problem.nControl;
nSegment = transcribe.param.nSegment;
nEven = transcribe.param.nEven;

x0.lb = problem.bnds.path0.lb;
x0.ub = problem.bnds.path0.ub;
xf.lb = problem.bnds.pathf.lb;
xf.ub = problem.bnds.pathf.ub;

F_lb = [ zeros( nSegment*nEven*nState , 1 ); zeros( nSegment , 1 ); x0.lb; xf.lb ];
F_ub = [ zeros( nSegment*nEven*nState , 1 ); zeros( nSegment , 1 ); x0.ub; xf.ub ];

end