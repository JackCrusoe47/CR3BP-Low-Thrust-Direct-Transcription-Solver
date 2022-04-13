function [solution,problem] = feasibleTranscription(problem,solver_options)
% =======================================================================
%       Direct Transcription/Collocation Feasible Trajectory Solver
%          Levenberg-Marquardt iterative least square algorithm
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 12-04-2022
%
% Format : [solution,problem] = ...
%              feasibleTranscription(problem,solver_options) 
%
% ***********************************************************************
%                             DESCRIPTION
% ***********************************************************************
%
%       A NLP solver to solve CR3BP+low-thrust trajectory through direct
%   transcription. It finds a feasible trajectory in proximity to a user
%   defined initial guess with optional boundary equality constraints.
%   Internally, the solver uses the Levenberg-Marquardt damped least
%   squares method with MATLAB fsolve solver.
% 
%   Solves the following trajectory optimal control problem:
%
%       given inital traj.  : xGuess, uGuess, tGuess
%
%       with dynamics       : dynamics(t,x,u,param)
%
%          & path const.    : u(1,t) <= maxThrust
%                             u(2,t)^2+u(3,t)^2+u(4,t)^2 = 1
%
%          & bnds const.    : x(t0) = x0    (OPTIONAL)
%                             x(tf) = xf    (OPTIONAL)
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
%       .path0.ub   : initial boundary equality constriants ( nState x 1 )
%       .pathf.ub   : final boundary equality costraints ( nState x 1 )
%
%   Set values for boundary equality constraints. For any unconstrainted 
%   states, set them to inf. For example, [5;2;3;inf;inf;inf;24], sets 
%   equality constraints for position as [5;2;3] and spacecraft mass as 24.
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
% -----------------------------------------------------------------------
%
% solver_options    : solver options struct
%
%   .maxIter        : maximum solver iterations (optional)
%   .tolFunc        : function tolerance (optional)
%   .tolStep        : step tolerance (optional)
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
%   .info           : output information struct from fsolve
%   .exitflag       : output exitflag from fsolve
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
%   .flag           : updated flag values
%       .path0      : initial state constraints flags ( nState x 1 )
%       .pathf      : final state constraints flags ( nState x 1 )
%
% ***********************************************************************
%
% ***********************************************************************
%                            CHANGE LOG
% ***********************************************************************
% 12-04-2022 : Code Created
% ***********************************************************************

%% feasible transcription code

% -- setup solver options

% - default fsolve options
solver_options.fsolve = optimoptions('fsolve','Display','iter',...
    'Algorithm','levenberg-marquardt',...
    'SpecifyObjectiveGradient',true,...
    'FunctionTolerance',1e-10,...
    'StepTolerance',1e-12,...
    'MaxIterations',50);

% - update max iterations
if isfield(solver_options,"maxIter")
    solver_options.fsolve.MaxIterations = solver_options.maxIter;
end

% - update function tolerance
if isfield(solver_options,"tolFunc")
    solver_options.fsolve.FunctionTolerance = solver_options.tolFunc;
    solver_options.fsolve.StepTolerance = solver_options.tolFunc * 1e-2; % two order lower than func. tol
end

% - update step tolerance (override)
if isfield(solver_options,"tolStep")
    solver_options.fsolve.StepTolerance = solver_options.tolStep;
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

% -- add boundary limits if absent
if ~isfield(problem.bnds,"path0")
    problem.bnds.path0.ub = +inf(problem.nState,1);
end
if ~isfield(problem.bnds,"pathf")
    problem.bnds.pathf.ub = +inf(problem.nState,1);
end

% -- check for path equality bounds 
problem.flag.path0 = zeros( problem.nState, 1 );
problem.flag.pathf = zeros( problem.nState, 1 );
problem.flag.path0( problem.bnds.path0.ub ~= +inf ) = 1;
problem.flag.pathf( problem.bnds.pathf.ub ~= +inf ) = 1;


% -- compute initial transcription from problem data
transcribe = initTrans(transcribe,problem);


% -- run IPOPT solver
solution = solveFeasible(transcribe,problem,solver_options);


% -- interpolate solution
solution = interpTraj(solution,transcribe,problem,solver_options);


% -- add transcription parameters
solution.param = transcribe.param;


end