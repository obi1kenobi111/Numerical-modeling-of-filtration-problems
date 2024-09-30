# Numerical-modeling-of-filtration-problems
This repository contains solutions to various problems of water and oil filtration in porous media on different branches.

All differential equations below are solved on uniform rectangular grids by the finite difference method with implicit time approximation. All equations except Task3 are linear and for them it is easy to obtain an asymmetric matrix (a system of linear equations). For a nonlinear equation in Task3 I use <a href="http://algowiki-project.org/ru/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%9D%D1%8C%D1%8E%D1%82%D0%BE%D0%BD%D0%B0_%D0%B4%D0%BB%D1%8F_%D1%81%D0%B8%D1%81%D1%82%D0%B5%D0%BC_%D0%BD%D0%B5%D0%BB%D0%B8%D0%BD%D0%B5%D0%B9%D0%BD%D1%8B%D1%85_%D1%83%D1%80%D0%B0%D0%B2%D0%BD%D0%B5%D0%BD%D0%B8%D0%B9" target="_blank">Newton's</a> method to obtain a system of linear algebraic equations (also an asymmetric matrix).
To solve the obtained asymmetric matrices I use the <a href="https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method" target="_blank">BiCGStab</a> method with the preconditioner <a href="https://en.wikipedia.org/wiki/Incomplete_LU_factorization" target="_blank">ILU</a>. The implementation of this matrix solving method and others can be found at 
<a href="https://github.com/kirill-terekhov/mipt-solvers?ysclid=m1p91idvpq689260986" target="_blank">Solvers</a> 

P.S: This repository contains my solutions to problems in the course "Numerical Methods in Oil and Gas Engineering" at Sirius University.

Examples of implict finite differences I used for linearization of differential equations:

<h3 align="center">${\left( \frac{d^2 f}{d x^2} \right)}^{n + 1} _ {ij}  \approx \frac{f^{n + 1} _{i-1,j} - 2 f^{n + 1} _{ij} + f^{n + 1} _{i + 1, j}}{\Delta x^2} $
  
<h3 align="center">${\left( \frac{df}{dt} \right)}^{n + 1} _ {ij} \approx \frac{f^{n + 1} _{i,j} - f^{n} _{ij}}{\Delta t} $
  
And <a href="https://en.wikipedia.org/wiki/Upwind_scheme" target="_blank">approximation</a> against the flow:

<h3 align="center">${\left( a \frac{df}{dx} \right)}^{n + 1} _ {ij} \approx max(a, 0) \frac{f^{n + 1} _{i,j} - f^{n + 1} _{i-1, j}}{\Delta x} + min(a, 0) \frac{f^{n + 1} _{i + 1,j} - f^{n + 1} _{i, j}}{\Delta x}$
    
Below are the formulations of these problems:

## branch Task1: 2D equation of single-phase (water) filtration in a square area
We consider a 2D equation of single-phase filtration in a square with a side of 1 (Ω=(0.1)×(0.1)) with a number of simplifying assumptions: incompressibility, neglecting the contribution of gravity, isotropic permeability tensor $K=1$. Then the equation will take a stationary form:

<h2 align="center">$\frac{\partial^2 p}{\partial x^2} + \frac{\partial^2 p}{\partial y^2} = q(x, y)$

I use Neumann boundary conditions everywhere: $\frac{\partial p}{\partial \mathbf{n}} = 0$ (no leakage)

The source member $q(x, y)$ models two wells: at point (0.25, 0.25) – injection well, pbh=1, at point (0.75, 0.75) – producing well, pbh=-1. To model wells I use the Peaceman formula:

<h2 align="center">$q(x, y) = \frac{2 \pi (p_{bh} - p)}{\ln(re/rw)}$

Discretization of second derivatives by the finite difference method is implicit (i.e. at time step $n + 1$):

<h2 align="center">${\left( \frac{d^2 f}{d x^2} \right)} _ {ij} \approx \frac{f_{i-1,j} - 2 f_{ij} + f{i + 1, j}}{\Delta x^2} $

So we get a system of linear equations.

## branch Task2: Non-stationary 2D water pressure equation for the single-phase filtration model in a square area

Water pressure equation:

<h2 align="center">$\frac{c_t \phi}{b_w}  \frac{\partial p}{\partial t}$ = $\frac{k}{\mu_w b_w} ( \frac{\partial^2 p}{\partial x^2} + \frac{\partial^2 p}{\partial y^2}) + q$

where

$b_w$ = 1 - incompressibility of a liquid

$\phi$ - porosity of the medium in which water is distributed

$k$ - permeability of the medium

$\mu_w$ - dynamic viscosity of water

The region contains one injection well at point $x = 5000 ft$, $y = 5000 ft$, fluid flow $q = 10^{−5} day^{−1}$ and one production well  at point $x = 9000 ft$, $y = 9000 ft$ with $p_bh = 800 psi$, $q = 2πln(re/rw)(p_{bh} − p)$
Initial pressure $p(t = 0) = 1000 psi$. The boundary conditions are the no-flow Neumann condition  $\frac{\partial p}{\partial \boldsymbol{n}} = 0$ on all boundaries except the right one $(x = L)$, where the Dirichlet condition is set $p = 2000 psi$

In this problem I also use implicit time discretization by the finite difference method

## branch Task3: Numerical and analytical solution of the Buckley-Leverett equation in 1D area

We will model the displacement of oil by water in a one-dimensional region by solving the equations to solve the Buckley-Leverett equation:

<h2 align="center">$\phi \frac{\partial S_w}{\partial t} = -\frac{q}{A} \frac{\partial f_w}{\partial x}$

where
$S_w, S_o$ - water saturation, oil saturation. $S_w$ + $S_o$ = 1 is performed for each point in the region

$f_w$ - water mobility, a nonlinear function of $S_w$ Thus, we get that our equations are nonlinear and you can’t just compose and solve the matrix of unknowns :That's why I implemented the solution of the nonlinear system using <a href="http://algowiki-project.org/ru/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%9D%D1%8C%D1%8E%D1%82%D0%BE%D0%BD%D0%B0_%D0%B4%D0%BB%D1%8F_%D1%81%D0%B8%D1%81%D1%82%D0%B5%D0%BC_%D0%BD%D0%B5%D0%BB%D0%B8%D0%BD%D0%B5%D0%B9%D0%BD%D1%8B%D1%85_%D1%83%D1%80%D0%B0%D0%B2%D0%BD%D0%B5%D0%BD%D0%B8%D0%B9" target="_blank">Newton's</a> method in the code.

$q = 1$ - Water flow at the left border of the region. So I will use the <a href="https://en.wikipedia.org/wiki/Upwind_scheme" target="_blank">scheme</a> with the difference against the fluxes for discretization of the spatial derivative, considering the flux $\frac{q}{A}$ to be positive.

The initial (at moment $t = 0$) saturation of water $S_w$ in the entire region is equal to its residual saturation $S_w = S_{wr}, S_{wr} = 0.2 = Const$. That is, the saturation of the environment with oil at $t = 0$ will be 0.8. Then water begins to flow into the region from the left border and displace the oil.

So, the equations are nonlinear because mobility $f_w$ depends nonlinearly on $S_w$:

$f_w$ = $\frac{m_w}{m_w + m_o}$ where $m_w = k_{rw}/\mu_w$ and $m_o = k_{ro}/\mu_o$ - local water and oil mobilities

The permeabilities of the phases ($k_{rw}$ and $k_{ro}$) are calculated using the Brooks-Corey model and are determined by the formula:

$k_{r \alpha}$ = $k_{r \alpha}^* S_{\alpha e}^{N_{\alpha}}$ $where$ $\alpha$ is phase $(w - water, or o - oil)$, $S_{\alpha e}$ = $\frac{S_\alpha - S_{\alpha r}}{1 - S_{or} - S_{wr}}$ - effective saturation

Parameters for the Brooks-Corey model in this task:

$k_{ro}^* = 1$

$k_{rw}^* = 0.6$

$N_o = N_w = 2$

$S_{wr} = 0.2, S_{or} = 0.15$

I also found an $analytical$  $solution$ for this problem, using the results from the <a href="https://www.researchgate.net/publication/309014212_Analysis_of_the_Buckley-Leverett_Solution_and_Comparison_with_Numerical_Simulation?enrichId=rgreq-0f016f82784ef3d9ae202b3bd53d4e8a-XXX&enrichSource=Y292ZXJQYWdlOzMwOTAxNDIxMjtBUzo0MTYzNDg0OTM2MzE0ODlAMTQ3NjI3NjYyNDE1MA%3D%3D&el=1_x_3&_esc=publicationCoverPdf" target="_blank">article</a>, which provides the dependence of the coordinates $x_D$ of a point with a fixed saturation $S_{wD}$:

<h2 align="center">$x_D = t_D \frac{d f_w}{d S_w} | _ {S_{wD}}  (*)$, $where$ $t_D$ = $\frac{q t}{\phi A}$

Obviously, the analytical solution for $S_w$ will look like this:

$S_w$ from (*) (we iterate over the values $S_w$ $\in$ $[1.0, S_{wf}]$), when $x <= x_f$ and $S_w = S_{wr}$, when $x > x_f$, where $x_f$ is wave front coorinate

To find $x_f$,I solved the equation $\frac{d f_w}{d S_w} | _ {S_{wf}}$ = $\frac{f_w (S_{wf})}{S_{wf} - S_{wr}}$ to find $S_{wf}$ - saturation of water on top of wave front and substitute it in (*)

## branch Task4: Single-phase flow in a 1D fractured porous formation
There is a crack running along the entire area. Single-phase flow in a fractured medium is described using the dual porosity/permeability model:

<h2 align="center">$\phi_{m} c_t \frac{\partial p_m}{\partial t}$ = $\frac{k_m}{\mu} \frac{\partial^2 p_m}{\partial x^2}$ + $\frac{k_m \lambda}{\mu} (p_f - p_m)$

<h2 align="center">$\phi_{f} c_t \frac{\partial p_f}{\partial t}$ = $\frac{k_f}{\mu} \frac{\partial^2 p_f}{\partial x^2}$ + $\frac{k_f \lambda}{\mu} (p_m - p_f)$

where $p_f$ - pressure in crack, $p_m$ - pressure in the area outside the crack

Thus, we have a system of two differential equations with unknowns $p_m$ and $p_f$

Boundry condition (BC) on the left boundary is $p(0, t) = 500$. BC on the right boundary – no leakage. Initial condition is $p(x, 0) = 1000$

# Output of results
The programs output results in <a href="https://www.paraview.org/" target="_blank">Paraview</a> (.vtk data files) and Excel (.csv data files) formats. The Paraview program is free and very quickly and easily installed in a couple of minutes. After installing paraview, select $File \to Open \to...$ to view the results.
