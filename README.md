# CR3BP Low-Thrust Direct Transcription Solver

A direct transcription solver (feasible estimator and optimizer) for low-thrust trajectories in circular restricted three-body problem.


## Features

* Direct transcription solver for CR3BP+low-thrust trajectories
* LGL node placement with user defined order
* Two different solvers for different stages of trajectory design
* Uses analytical jacobians for fast computation.
* Uses sparse matrices to handle very large problems.
* Feasible trajectory solver
	* Estimate closest feasible trajectory from an initial guess 
	* Levenberg-Marquardt with MATLAB fsolve solver
* Optimizer
	* Optimize input trajectory with various objectives (thrust, propellant and time optimal)
	* Allows robust feasible trajectory estimation with dummy zero objective
	* Interior-Point with IPOPT optimizer


## Requirements

* [**LGL nodes and weights**][1] by Greg von Winckel
* [**mexIPOPT**][2] IPOPT optimizer port for MATLAB by Enrico Bertolazzi

	[1]: https://www.mathworks.com/matlabcentral/fileexchange/4775-legende-gauss-lobatto-nodes-and-weights
	[2]: https://github.com/ebertolazzi/mexIPOPT

## Installation

* Download CR3BP Low-Thrust Direct Transcription Solver
* Download and install IPOPT optimizer
* Download LGL nodes-and-weights and copy to **external\\** folder

## Limitation

* No straight forward method exists for adding user-defined constraints, except what is already implemented. Currently, the code only supports the two-point boundary-value problem for CR3BP+low-thrust. Adding other constraints (for example periodicity) requires cubersome modification of the code. 
* Some parts of the code are developed with modularity and abstraction in mind. For example, the transcription defect constraints and jacobain can be computed for any user defined dynamics provided the analytical jacobian wrt. state, control and time are available. But most other sections are tailored for CR3BP+low-thrust problem.
* Through manual mesh refinement can be achived by chaining multiple solver instances with varying segment and order values externally, no inbuilt mesh refinement routines currently are implemented.
* For use with large problems, the code uses sparse matrices for the full jacobian of the NLP. But sparse matrix creation currently is not fully optimized and code readablity is terrible in its implentation.
