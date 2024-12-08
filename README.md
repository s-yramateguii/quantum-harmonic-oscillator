# Quantum Harmonic Oscillator

This repository contains Python code that solves for the eigenfunctions and eigenvalues of the quantum harmonic oscillator using various numerical methods. The code implements both linear and nonlinear solvers, and compares methods like the shooting method and finite difference methods. Additionally, the performance of different ODE solvers is analyzed.

## Table of Contents

- [Overview](#overview)
- [Mathematical Formulation](#mathematical-formulation)
- [Methods Used](#methods-used)
  - [Shooting Method](#shooting-method)
  - [Finite Difference Method](#finite-difference-method)
  - [Nonlinear Solvers](#nonlinear-solvers)

## Overview

This project solves the time-independent Schrödinger equation for a quantum harmonic oscillator in both linear and nonlinear settings. It uses different numerical methods to find eigenvalues and eigenfunctions. The eigenfunctions are normalized and the errors in the computed eigenvalues and eigenfunctions are compared to exact results.

The problems tackled include:

- **Part (a):** Solving for first five eigenfunctions and eigenvalues using the shooting method.
- **Part (b):** Solving for first five eigenfunctions and eigenvalues using the finite difference method and comparison with the shooting method.
- **Part (c):** Solving a nonlinear Schrödinger equation for first two eigenfunctions and eigenvalues using shooting method.
- **Part (d):** Analyzing the performance of various ODE solvers with different tolerance settings.
- **Part (e):** Comparing the errors in eigenfunctions and eigenvalues obtained from both methods with exact solutions.

## Mathematical Formulation

The primary equation solved in this project is the Schrödinger equation:

$$
\frac{d^2 \phi_n}{dx^2} + [Kx^2-\epsilon_n] \phi_n = 0
$$

Where:
- $$\phi(x)$$ is the wavefunction.
- $$\epsilon_n$$ is the quantum energy.
- $$K$$ is the spring constant (or stiffness) for the harmonic oscillator ($$K=1$$ for all the parts).

The eigenfunctions $$\phi(x)$$ are solutions to this equation for discrete energy eigenvalues, and they can be expressed in terms of Hermite polynomials:

$$
\phi(x) = N_n H_n(x) e^{-x^2 / 2}
$$

Where:
- $$N_n$$ is the normalization constant.
- $$H_n(x)$$ is the $n$-th Hermite polynomial.
 
The nonlinear equation is also used:

$$
\frac{d^2 \phi_n}{dx^2} - [\gamma|\phi_n|^2+Kx^2-\epsilon_n] \phi_n = 0
$$

Where:
- $$\gamma$$ is the probability density
## Methods Used

### Shooting Method

In this approach, the eigenvalue $$\epsilon_n$$ is adjusted iteratively using a boundary condition check at the boundaries of the domain. The shooting method involves solving the differential equation for different guesses of $$\epsilon_n$$ and adjusting $$\epsilon_n$$ until the solution satisfies the boundary conditions. This is done for both linear and nonlinear Schrödinger equations.

### Finite Difference Method

A second-order finite difference method is used to discretize the differential equation, which results in a matrix eigenvalue problem. The matrix is then diagonalized to obtain the eigenvalues and eigenfunctions.

### Nonlinear Solvers

For the nonlinear Schrödinger equation, a similar shooting method is applied, but with the nonlinear term $$\gamma |\phi(x)|^2 \phi(x)$$ included in the differential equation. The energy eigenvalue $$\epsilon_n$$ is adjusted in a similar fashion, but the solution now includes nonlinearity. $$\gamma=\pm0.05$$ is used.
