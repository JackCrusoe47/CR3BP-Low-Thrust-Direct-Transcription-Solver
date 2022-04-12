function transcribe = initTrans(transcribe,problem)
% -- initial trajectory transcription from problem data

% - segment size
Deltat = ( problem.traj.data.t(end)-problem.traj.data.t(1) )/transcribe.param.nSegment;
transcribe.traj.Deltat = Deltat;

% - vector transcription time nodes
t = 0;
for n = 1:transcribe.param.nSegment
    for i = 2:transcribe.param.nOdd
        t = [t, (n-1)*Deltat + (transcribe.param.tauOdd(i)+1)*Deltat/2];
    end
end

% - interpolate from problem data
transcribe.traj.nodes.t = t;
transcribe.traj.nodes.x = problem.traj.interp.x(t);
transcribe.traj.nodes.u = problem.traj.interp.u(t);

end