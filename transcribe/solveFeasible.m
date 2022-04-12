function solution = solveFeasible(transcribe,problem,solver_options)

nState = problem.nState;
nControl = problem.nControl;
nSegment = transcribe.param.nSegment;
nNodes = length(transcribe.traj.nodes.t);
nOdd = transcribe.param.nOdd;
timeFixed = problem.flag.timeFixed;

% -- initial solution for optimization variable
z = real(optimVar_feasible(transcribe,problem));

% Set up the auxiliary data.
auxdata.transcribe = transcribe;
auxdata.problem = problem;

% -- setup fsolve problem
data.objective = @(z)solveFunc(z,auxdata);
data.x0 = z;
data.solver = 'fsolve';
data.options = solver_options.fsolve;

[zsol,fval,exitflag,output]  = fsolve(data);

% -- solution
solution.info = output;
solution.fval = fval;
solution.exitflag = exitflag;
solution.zVar.initial = z;
solution.zVar.solved = zsol;

if ~timeFixed
    lenz = (nSegment*(nOdd-1)+1)*nState + nSegment*nControl + 1;
else
    lenz = (nSegment*(nOdd-1)+1)*nState + nSegment*nControl;
end

zMain = zsol(1:lenz);

if ~timeFixed
    solution.traj.nodes.x = reshape( zMain(1:nNodes*nState), nState, nNodes);
    solution.traj.nodes.u = reshape( zMain(nNodes*nState+1:end-1), nControl, nSegment);
    solution.traj.nodes.t = (transcribe.traj.nodes.t) .* zMain(end)/transcribe.traj.nodes.t(end);
else
    solution.traj.nodes.x = reshape( zMain(1:nNodes*nState), nState, nNodes);
    solution.traj.nodes.u = reshape( zMain(nNodes*nState+1:end), nControl, nSegment);
    solution.traj.nodes.t = transcribe.traj.nodes.t(end);
end

    function [F,DF] = solveFunc(z,auxdata)
        
        % - get constraints vector
        F = constFunc_feasible(z,auxdata);
        
        % - get jacobian matrix
        DF = jacFunc_feasible(z,auxdata);
        
    end

end