# Verification-of-a-counterexample-given-to-junction-points-optimality-conditions-in-singular-control

This repository contains the MATLAB script used as a companion to the research paper **"On Mcdanell's Conjecture"**. 

## Description

The code verifies that the counterexample given in the paper **"Junction Points in Singular Optimal Control"** is indeed a valid counterexample using MATLAB Symbolic Toolbox. It calculates the values of:

- \(\alpha^{(m)}(0)\) and \(\beta^{(m)}(0)\) for \(m = 0, 1, 2\)
- \(u_s^{(m)}(0)\) for \(m = 0, 1, 2, 3\)

and compares them to the values provided in the paper.

## Usage

The code can be executed directly, and it will display the results by printing the steps and underlying logic.

### Note on Performance:
- Calculating the derivatives using MATLAB Symbolic Toolbox can take some time.
- For faster verification, precomputed results are provided in a saved workspace file. Users can load this file instead of running the computations:
  ```matlab
  load('results.mat');
