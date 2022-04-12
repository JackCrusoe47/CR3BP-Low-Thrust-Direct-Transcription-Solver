function grad = gradFunc_timeOptimal(z,transcribe,problem)

% % - extract parameters
% nState = problem.nState;
% nSegment = transcribe.param.nSegment;
% nOdd = transcribe.param.nOdd;
timeFixed = problem.flag.timeFixed;

% - gradient of objective
grad = zeros(size(z));
if ~timeFixed
    grad(end) = 1;
end

end