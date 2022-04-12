function [tau,W,P]=nodes_LGL(N)
% function wrapper for LGL nodes algorithm by Greg von winckel
%
% Ref: Greg von Winckel (2022). Legende-Gauss-Lobatto nodes and weights 
%      (https://www.mathworks.com/matlabcentral/fileexchange/4775-legende
%      -gauss-lobatto-nodes-and-weights), MATLAB Central File Exchange. 
%      Retrieved April 11, 2022.
%
% Download the file or uncomment if original file not found. All credits of
% LGL node computation code is to Greg von winckel.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  COMMENT ONLY IF ORIGIANAL FILE NOT AWAILABLE IN MATLAB FILE EXCHANGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tau,W,P]=lglnodes(N-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  UNCOMMENT ONLY IF ORIGIANAL FILE NOT AWAILABLE IN MATLAB FILE EXCHANGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % order correction.
% N = N-1;
% 
% % Truncation + 1
% N1=N+1;
% % Use the Chebyshev-Gauss-Lobatto nodes as the first guess
% tau=cos(pi*(0:N)/N)';
% % The Legendre Vandermonde Matrix
% P=zeros(N1,N1);
% % Compute P_(N) using the recursion relation
% % Compute its first and second derivatives and 
% % update x using the Newton-Raphson method.
% xold=2;
% while max(abs(tau-xold))>eps
%     xold=tau;
%         
%     P(:,1)=1;    P(:,2)=tau;
%     
%     for k=2:N
%         P(:,k+1)=( (2*k-1)*tau.*P(:,k)-(k-1)*P(:,k-1) )/k;
%     end
%      
%     tau=xold-( tau.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
%              
% end
% W=2./(N*N1*P(:,N1).^2);
% 

% flipping vectors to be inline with transcription routine
tau = flip(tau);
W = flip(W);