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
%   A feasible trajectory solver for CR3BP+low-thrust problem with direct
% transcription. It estimates a feasible / quasi-feasible trajectory in 
% proximity of an initial guess solution provided by the user through
% through iterative correction with a Newton-Raphson like scheme. The
% solver can impose optional boundary state constraints along with 
% transcription defects. For robustness, instead of a perfect Newton - 
% Raphson apporoach, the method used is the Levenberg-Marquardt damped 
% least square method. 
%
%   For accurate transcription, instead of unifrom nodes, the nodes are
% distributed using a Legandre-Gauss-Lobotto (LGL) node placement. To
% handle large problems, the jacobian matrices of constraints w.r.t.
% to the variables are implemented with sparse matrix, saving valuble
% memorey. For efficiency, analytical jacobians are computed instead of
% expensive finite deference methods. For simplicity and being inline with
% typical low-thrust maneuvers, the controls are kept constant within a
% segment.
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
% For current CR3BP+LowThrust, set these functions to dynamics_CR3BP_Thrust
% and jacobian_CR3BP_Thrust respectively.
%
%   .bnds           : problem bounds struct (OPTIONAL)
%       .path0.ub   : initial boundary equality constriants ( nState x 1 )
%       .pathf.ub   : final boundary equality costraints ( nState x 1 )
%
% Set values for boundary equality constraints. For any unconstrainted 
% states, set them to inf. For example, [5;2;3;inf;inf;inf;24], sets 
% equality constraints for position as [5;2;3] and spacecraft mass as 24.
%                      
%   .flag           : additional problem flags
%       .timeFixed  : true if problem is time-fixed (not recommended)
%
% The time fixed option is not recommended. For most problems, it prevents 
% the solver from converging, either due to overconstraining or error in 
% jacobian.s
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
% Content same as input problem struct, along with updates listed below:
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