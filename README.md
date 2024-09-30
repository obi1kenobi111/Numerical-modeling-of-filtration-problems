# Numerical-modeling-of-filtration-problems
This repository contains solutions to various problems of water and oil filtration in porous media on different branches. Below are the formulations of these problems:

## branch Task1: 2D equation of single-phase filtration in a square area
We consider a 2D equation of single-phase filtration in a square with a side of 1 (Ω=(0.1)×(0.1)) with a number of simplifying assumptions: incompressibility, neglecting the contribution of gravity, isotropic permeability tensor $K=1$. Then the equation will take a stationary form:

$\frac{\partial^2 p}{\partial x^2} + \frac{\partial^2 p}{\partial y^2} = q(x, y)$

I use Neumann boundary conditions everywhere: $\frac{\partial p}{\partial \boldsymbol{n}} = 0$

The source member $q(x, y)$ models two wells: at point (0.25, 0.25) – injection well, pbh=1, at point (0.75, 0.75) – producing well, pbh=-1. To model wells I use the Peaceman formula:

$q(x, y) = \frac{2 \pi \(p_{bh} - p)}{\ln(re/rw)}$

Discretization of second derivatives by the finite difference method is implicit (i.e. at time step $ n + 1 $):

${\left( \frac{d^2 f}{d x^2} \right)} _ {ij} \approx \frac{f_{i-1,j} - 2 f_{ij} + f{i + 1, j}}{\Delta x^2} $

So we get a system of linear equations.

## branch Task2: 2D equation of single-phase filtration in a square area

$\frac{c_t \phi}{b_w}  \frac{\partial p}{\partial t}$ = $\frac{k}{\mu_w b_w} ( \frac{\partial^2 p}{\partial x^2} + \frac{\partial^2 p}{\partial y^2}) + q$
