function [solution,problem] = optimalTranscribtion(problem,solver_options)
% =======================================================================
%         Direct Transcription/Collocation Trajectory Optimizer
%               IPOPT interior point optimization method
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 12-04-2022
%
% Format : [solution,problem] = ...
%              optimalTranscribtion(problem,solver_options) 
%
% ***********************************************************************
%                             DESCRIPTION
% ***********************************************************************
%
%       A NLP optimizer to solve CR3BP+low-thrust trajectory through direct
%   transcription. It finds an optimized trajectory in proximity to a user
%   defined initial guess with given a objective and inequality constraints
%   on boundary & path parameters. The optimization is performed with an
%   interior point optimizer IPOPT.
% 
%   Optimizes the following trajectory optimal control problem:
%
%       given objective     : min( J )
%
%          & inital traj.   : xGuess, uGuess, tGuess
%
%       with dynamics       : dynamics(t,x,u,param)
%
%          & path const.    : xlb <= x(t) <= xub    (OPTIONAL)
%                             uln <= u(t) <= uub    (OPTIONAL)
%                             u(1,t) <= maxThrust
%                             u(2,t)^2+u(3,t)^2+u(4,t)^2 = 1
%
%          & bnds const.    : x(t0) = x0    (OPTIONAL)
%                             x(tf) = xf    (OPTIONAL)
%                             tlb <= tf <= tub (OPTIONAL) 
%
%       Using the 'feasible' dummy objective, the optimizer can solve NLP 
%   to identify a feasible trajectory much more robustly than the feasible
%   transcription solver.
%
%       Converts continous optimal control problem to a NLP with direct
%   transcription method. The trajectory is descritized in to segments
%   and special nodes in each segment are used to create an interpolation
%   function and defects in collocation nodes are minimized. For improved
%   accuracy, Legandre-Gauss-Lobotto node placement strategy is used to
%   define the node location. The user defines the number of segments and 
%   order of interpolation for the transcription. The program uses 
%   analytical jacobian of the dynamics and uses sparse matrices for
%   computational and memory efficency.
%
% ***********************************************************************
%
% ***********************************************************************
%                               INPUTS
% ***********************************************************************
% 
% problem           : input problem struct
%
%   .param          : problem parameters substruct
%    	.mu         : CR3BP characteristic parameter [nd]
%    	.tstr       : CR3BP characteristic time [s]
%    	.lstr       : CR3BP characteristic length [km]
%    	.Isp(Thr)   : Specific Impulse with thrust function (Thr) [s]
%    	.dIsp(Thr)  : Specific Impulse derivative with Thrust (Thr) [s/N]
%       .maxThrust  : Maximum thrust level [N]
%       .maxMass    : Maximum spacecraft mass [kg]
%
%   .traj           : input trajectory struct
%       .data       : trajectory raw data struct
%           .t      : time vector ( 1 x nPoints )
%           .x      : state vector ( nState x nPoints )
%           .u      : control vector (nControl x nPoints )
%
%   .func           : problem function handles
%       .dynamics   : dynamics function handle @(t,x,u)dynProb(t,x,u,param)
%       .jacobian   : jacobian function handle @(t,x,u)jacProb(t,x,u,param)
%
%   For current CR3BP+LowThrust, set these functions to 
%   dynamics_CR3BP_Thrust and jacobian_CR3BP_Thrust respectively.
%
%   .bnds           : problem bounds struct (OPTIONAL)
%       .state      : state bounds (OPTIONAL)
%           .lb     : lower bound ( nState x 1 )
%       	.ub     : upper bound ( nState x 1 )
%       .control    : control bounds (OPTIONAL)
%           .lb     : lower bound ( nControl x 1 )
%       	.ub     : upper bound ( nControl x 1 )
%       .time       : time bounds (OPTIONAL)
%           .lb     : lower bound ( 1 x 1 )
%           .ub     : upper bound ( 1 x 1 )
%       .path0      : initial state bounds (OPTIONAL)
%       	.lb     : lower bound ( nState x 1 )
%       	.ub     : upper bound ( nState x 1 )
%       .pathf      : final state bounds (OPTIONAL)
%       	.lb     : lower bound ( nState x 1 )
%       	.ub     : upper bound ( nState x 1 )
%
%   Recommended to leave all, but the most necessary option, to their 
%   defaults. The rate of convergence can be reduced considerably with over
%   constrained or redundant constraints.
%                      
%   .flag           : additional problem flags
%       .timeFixed  : true if problem is time-fixed (not recommended)
%
%   The time fixed option is not recommended. For most problems, the option
%   prevents the solver from converging, either due to overconstraining or 
%   error in jacobian.
%
%   .nState         : number of states (7 for CR3BP+Low Thrust problem)
%   .nControl       : number of controls (4 for CR3BP+Low Thrust problem)
%
%   .objType        : objective of optimization ( string )
%
%                     'thrustOptimal'   : minimize( sum(u(t)*u(t)) )
%                     'propOptimal'     : maximize( m(tf) )
%                     'timeOptimal'     : minimize( tf )
%                     'stateTarget'     : minimize( x(tf) - xTarg )
%                     'massTarget'      : minimize( m(tf) - mTarg )
%                     'feasible'        : minimize( 0 ) (DEFAULT)
%
% -----------------------------------------------------------------------
%
% solver_options    : solver options struct
%
%   .maxIter        : maximum solver iterations (optional)
%   .tolFunc        : function tolerance (optional)
%
%   .nSegment       : number of transcription segments
%   .nOrder         : order of transcription (odd number)
%
%   .dtInterp       : output interpolation step size
%
% ***********************************************************************
%
% ***********************************************************************
%                              OUTPUTS
% ***********************************************************************
%
% solution          : output solution struct
%
%   .traj           : output trajectory struct
%
%       .data       : 
%           .t      : output interpolated time points ( 1 x nOutput )
%           .x      : output interpolated states ( nState x nOutput )
%           .u      : output interpolated control ( nControl x nOutput )
%
%       .nodes      : raw interpolation nodes used for transcription
%           .t      : time points at nodes ( 1 x nNodes )
%           .x      : states at nodes ( nState x nNodes )
%           .u      : controls at segments ( nControl x nSegment )
%
%       .interp     : interpolation functions using traj.data
%           .x(t)   : state interpolation function ( nState x nTime )
%           .u(t)   : control interpolation function ( nControl x nTime )
%
%   .zVar           : raw solver variables used for transcription
%       .initial    : initial values of variables ( nVar x 1 )
%       .solved     : final values of variables ( nVar x 1 )
%
%   .fval           : constraint vector at solution ( nConst x 1 )
%   .info           : output information struct from IPOPT
%   .exitflag       : output exit status from IPOPT
%   .objective      : objective value from IPOPT
%
% -----------------------------------------------------------------------
%
% problem           : updated problem struct for debuging
%
%   Content same as input problem struct, along with updates listed below:
%
%   .traj.interp    : interpolation functions for input traj.data
%           .x(t)   : state interpolation function ( nState x nTime )
%           .u(t)   : control interpolation function ( nControl x nTime )
%
% ***********************************************************************
%
% ***********************************************************************
%                            CHANGE LOG
% ***********************************************************************
% 12-04-2022 : Code Created
% ***********************************************************************

%% optimizer transcription code

% -- setup solver options

% - IPOPT basic options
solver_options.ipopt.jac_d_constant   = 'no';
solver_options.ipopt.hessian_constant = 'no';
solver_options.ipopt.mu_strategy      = 'adaptive';
solver_options.ipopt.max_iter         = 100;
solver_options.ipopt.tol              = 1e-10;

% - hessian computation options
solver_options.ipopt.hessian_approximation      = 'limited-memory';
solver_options.ipopt.limited_memory_update_type = 'bfgs';

% - update max iterations
if isfield(solver_options,"maxIter")
    solver_options.ipopt.max_iter = solver_options.maxIter;
end

% - update function tolerance
if isfield(solver_options,"tolFunc")
    solver_options.ipopt.tol = solver_options.tolFunc;
end


% -- compute transcription parameters
transcribe.param = getParamTranscribe_LGL(solver_options.nOrder);
transcribe.param.nSegment = solver_options.nSegment;

% -- add interpolation to problem initial solution
problem.traj.interp.x = @(z) interp1(problem.traj.data.t',problem.traj.data.x',z','pchip')';
problem.traj.interp.u = @(z) interp1(problem.traj.data.t',problem.traj.data.u',z','pchip')';

% -- add bnds substruct if absent
if ~isfield(problem,"bnds")
    problem.bnds = struct();
end

% -- add state limits if absent
if ~isfield(problem.bnds,"state")
    problem.bnds.state.lb = -inf(problem.nState,1);
    problem.bnds.state.ub = inf(problem.nState,1);
end

% -- add control limits if absent
if ~isfield(problem.bnds,"control")
    problem.bnds.control.lb = [0;-inf;-inf;-inf];
    problem.bnds.control.ub = [problem.param.maxThrust;inf;inf;inf];
end

% -- add time limits if absent
if ~isfield(problem.bnds,"time")
    problem.bnds.time.lb = 0;
    problem.bnds.time.ub = +inf;
end

% -- add boundary limits if absent
if ~isfield(problem.bnds,"path0")
    problem.bnds.path0.lb = [-inf(6,1);0];
    problem.bnds.path0.ub = [+inf(6,1);problem.param.maxMass];
end
if ~isfield(problem.bnds,"pathf")
    problem.bnds.pathf.lb = [-inf(6,1);0];
    problem.bnds.pathf.ub = [+inf(6,1);problem.param.maxMass];
end

% -- add dummy objective if absent
if ~isfield(problem,"objType")
    problem.objType = 'feasible';
end


% -- compute initial transcription from problem data
transcribe = initTrans(transcribe,problem);


% -- run IPOPT solver
solution = solveIPOPT(transcribe,problem,solver_options);


% -- interpolate solution
solution = interpTraj(solution,transcribe,problem,solver_options);


% -- add transcription parameters
solution.param = transcribe.param;


% -- add additonal solution outputs
solution.fval = solution.info.lambda;
solution.exitflag = solution.info.status;
solution.objective = solution.info.objective;

end