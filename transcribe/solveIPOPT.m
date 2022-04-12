function solution = solveIPOPT(transcribe,problem,solver_options)

nState = problem.nState;
nControl = problem.nControl;
nSegment = transcribe.param.nSegment;
nNodes = length(transcribe.traj.nodes.t);
timeFixed = problem.flag.timeFixed;

% -- initial solution for optimization variable
z = optimVar(transcribe,problem);

% -- optimization variable bounds
[z_lb,z_ub] = bndsVar(transcribe,problem);
options.lb = z_lb' ;    % Lower bound on the variables.
options.ub = z_ub' ;    % Upper bound on the variables.

% -- optimization constraint bounds
[F_lb,F_ub] = bndsConst(transcribe,problem);
options.cl = F_lb';     % Lower bound on the constraints.
options.cu = F_ub';     % Upper bound on the constraints.

% Set up the auxiliary data.
auxdata.transcribe = transcribe;
auxdata.problem = problem;
options.auxdata = auxdata ;

% Set the IPOPT options.
options.ipopt = solver_options.ipopt;

% The callback functions.
funcs.objective         = @objFunc;
funcs.constraints       = @constFunc;
funcs.gradient          = @gradFunc;
funcs.jacobian          = @jacFunc;
funcs.jacobianstructure = @strFunc;

[zopt, info] = ipopt_auxdata(z,funcs,options);

% -- solution
solution.info = info;
solution.zVar.initial = z;
solution.zVar.optimal = zopt;
if ~timeFixed
    solution.traj.nodes.x = reshape( zopt(1:nNodes*nState), nState, nNodes);
    solution.traj.nodes.u = reshape( zopt(nNodes*nState+1:end-1), nControl, nSegment);
    solution.traj.nodes.t = (transcribe.traj.nodes.t) .* zopt(end)/transcribe.traj.nodes.t(end);
else
    solution.traj.nodes.x = reshape( zopt(1:nNodes*nState), nState, nNodes);
    solution.traj.nodes.u = reshape( zopt(nNodes*nState+1:end), nControl, nSegment);
    solution.traj.nodes.t = transcribe.traj.nodes.t(end);
end


end