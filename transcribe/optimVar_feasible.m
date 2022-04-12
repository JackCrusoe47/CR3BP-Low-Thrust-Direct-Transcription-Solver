function z = optimVar_feasible(transcribe,problem)

% - extract parameters
nSegment = transcribe.param.nSegment;
nOdd = transcribe.param.nOdd;
maxThrust = problem.param.maxThrust;

% - normal transcription variables
z = optimVar(transcribe,problem);

% - control vector for segments
u = transcribe.traj.nodes.u(:, 1:nOdd-1:nSegment*(nOdd-1));

% - slack variables on thrust
zSlack = asin( ( u(1,:)/maxThrust ).^(1/2) );

z = [z; zSlack'];

end