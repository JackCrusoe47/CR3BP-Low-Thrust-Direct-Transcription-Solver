# CR3BP Low-Thrust Direct Transcription Solver

A direct transcription solver (feasible estimator and optimizer) for low-thrust trajectories in circular restricted three-body problem.


## Features

* Direct transcription solver for CR3BP+low-thrust trajectories
* LGL node placement with user defined order
* Two different solvers for different stages of trajectory design
* Uses analytical jacobians for fast computation.
* Uses sparse matrices to handle very large problems.
* Feasible trajectory solver (**feasibleTranscription()**)
	* Estimate closest feasible trajectory from an initial guess 
	* Levenberg-Marquardt with MATLAB fsolve solver
* Optimizer (**optimalTranscription()**)
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

## Usage / Manual

* Use ```help feasibleTranscription``` for brief information on feasible solver.
* Use ```help optimalTranscription``` for brief information on optimizer.
* Examples of solver usage can be found in **examples\\** folder.

## Limitation / Issues

* No direct access to LGL interpolation outside solver. 
* No easy way to add user-defined constraints (eg: periodicity)
* Optimization issues with some sparse matrix operations
* Time fixed flag is buggy
* No inequality constraints for feasible transcription solver
* Small part of the code is built with modularity in mind, others are tailored for CR3BP+low-thrust

## References

- "Low-Thrust Trajectory Design for Tours of the Martian Moons", Beom Park
- "A Collocation Approach for Computing Solar Sail Lunar Pole-Sitter Orbits", Martin T. Ozimek, Daniel J. Grebow, and Kathleen C. Howell
- "MColl - Monte Collocation Trajectory Design Tool", Daniel J. Grebow and Thomas A. Pavlak
- "Strategies for Low-Thrust Transfer Design Based on Direct Collocation Techniques", Robert E. Pritchett
- "Survey of Direct Transcription for Low-Thrust Space Trajectory Optimization with Applications", Francesco Topputo and Chen Zhang
- "Trajectory Design in the. Earth-Moon System and Lunar South Pole Coverage", Daniel J. Grebow
- "Design of low-thrust transfers. from an NRHO to low lunar orbits-applications for small spacecraft", Beom Park, Kathleen C. Howell, and David C. Folta
