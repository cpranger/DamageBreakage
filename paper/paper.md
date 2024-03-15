Numerical modeling of multi-state rock rheology
=====

Casper Pranger, Dave May, Alice Gabriel, et al.

Abstract
--------
We present a numerical scheme for efficiently modeling the rheological evolution of rocks in response to changes in temperature, damage degree, and grain size. The intended range of application includes evolution towards rheological instability and subsequent recovery in lab- or crustal scale domains of up to three spatial dimensions. To this end, elastodynamic effects are modeled as required by the rheological time scale. We achieve these aims at reasonable computational cost through the use of an additive implicit-explicit (IMEX) split Runge-Kutta scheme with error estimation and step size control, in which individual PDE terms are solved implicitly, and coupling terms are solved explicitly. This results in a method in which numerical stability is not limited by spatial resolution.


# Table of contents
- 1: Introduction
- 2: Problem statement [Merge into 1?]
- 3: Numerical methods
	- 3.1: Considerations for time discretization
	- 3.2: An additive Runge-Kutta method
	- 3.3 Stability and Accuracy
	- 3.4: Error estimation and time step control
	- 3.5: Decomposition of the governing equations
	- 3.6: Iterative non-linear solution methods
	- 3.7: Spatial discretization
- 4: Modeling results
	- 1D Thermo-elasticity
	- 1D Nonlinear wave propagation (cf. Zihua)
	- 2D Thermal runaway + grain size (cf. Thielmann)
	- 2D Damage + breakage (cf. Ben-Zion)
- 5: Discussion
- 6: Conclusion

-----------------------


1: Introduction
---------------

2: Problem statement
-----------------

<!-- Let $t \in \mathrm{T}$ denote a coordinate in a half-open interval $\mathrm{T} = [0,l_0) \subset \mathbb{R}$ of time with boundary $\partial \mathrm{T} = \{0\}$, and
let $x \in \Omega$ be a coordinate in a rectangular region $\Omega = [0,l_1] \otimes [0,l_2] \otimes [0,l_3] \subset \mathbb{R}^3$ of 3-space with boundary $\partial\Omega$. -->

We wish to examine the behavior of an elastic solid occupying the region $\mathrm{T} \otimes \Omega$ with its initial configuration specified on $\partial \mathrm{T}$, that is steadily loaded by displacement boundary conditions on $\partial\Omega_u \subseteq \partial\Omega$, and occasionally rapidly unloaded by transient anelastic processes, by transient material degradation, or by combinations thereof. To this end we consider Cauchy's linear momentum balance law

$$\tag{eq:2.1}
    d_t v = r^{-1} \nabla \cdot s + g,
$$

together with the mass balance law

$$\tag{eq:2.2}
    d_t r = -r \nabla \cdot v,
$$

in which $r = r(t,x(t))$ denotes the mass density at time $t$ at material coordinate $x(t)$, $v = v(t,x(t))$ likewise the material velocity, $s = s(t,x(t))$ Cauchy's stress tensor, and $g = g(t,x(t))$ represents the action of body forces, gravitational or otherwise. Body forces are ignored from here on ($g = 0$) at no loss of generality -- although the physical framework developed here ignores self-gravitation. At considerable loss of generality, but required for the sake of concise argument we ignore advection and distortion within the elastic body by setting $x(t) = x$, with $x$ the coordinate in the reference configuration. The total time derivative $d_t f(t, x^\ast(t))$ along the characteristic of flow $x(t)$ is consequently replaced by the partial time derivative $\partial_t f(t, x)$. The divergence $\nabla \cdot s$ of the Cauchy stress tensor $s$, which formally is to be evaluated in deformed coordinates $x(t)$, is instead evaluated in the reference coordinate system $x$. The density deviations due to volume change are assumed to have a negligible impact on the dynamic momentum balance, i.e. $d_t r = \partial_t r = 0$.

<!-- The displacement field $u = u(t,x)$ is governed by

$$\tag{eq:2.3}
    \frac{\tilde{d}u}{dt} = v \approx \frac{\partial u}{\partial t}.
$$ -->

We will revisit the implications of these choices in the context of the present methodological work in the discussion in Section [Discussion].

<!-- %Combining the mass balance law {eq:2.1} with {eq:2.3} and solving for $r = r(u)$ gives
% \begin{align} \label{eq:density}
	% r(u) = r_0 \exp \nabla \cdot (u - u_0) = r_0 \exp \nabla \cdot u,  \quad u(t = 0,x) &= u_0(x) = 0, \\
	% r(t=0,x) &= r_0(x). \nonumber
% \end{align}
% This solution eliminates differential equation {eq:2.1} and unknown $r$ from the problem. -->

In line with our assumption of a sufficiently distortion-free elastic medium is the use of the (symmetric) infinitesimal strain tensor $\varepsilon = \varepsilon(t, x)$, which is additively decomposed into an elastic part $e = e(t,x)$ and an anelastic part $\bar{e} = \bar{e}(t,x)$ and relates to the displacement $u = u(t, x)$ as

$$\tag{eq:2.4}
	\varepsilon = \nabla^\mathrm{s} u = \frac{1}{2} \left[ (\nabla u)^\mathrm{T} + (\nabla u) \right] = e + \bar{e}.
$$

<!-- % We choose to invert the relationship {eq:2.4} for the function
% \begin{align} \label{eq:inftalstrain2}
	% e_\mathrm{e} = e_\mathrm{e}(\nabla^\mathrm{s} u, \bar{e}) = \nabla^\mathrm{s} u - \bar{e}.
% \end{align} -->

An isotropic but otherwise quite general nonlinear elastic constitutive equation with structural stresses [e.g. Lyakhovksy et al., 2014b] is given by

$$\tag{eq:2.5}
	s(e, \alpha) = \lambda(x,\alpha, \jmath e) \,\delta \;\mathrm{tr}\; e + 2 \mu(x,\alpha, \jmath e) e - \vartheta(x) (\nabla \alpha) \otimes (\nabla \alpha).
$$

Here $\lambda$ and $\mu$ are the so-called Lam\'e parameters of an isotropic elastic solid, and $\vartheta$ an additional modulus associated with the structural stresses incurred by the damage field $\alpha = \alpha(t,x)$. The notation $\jmath e$ refers to the three scalar invariants $(J_1, J_2, J_3)$ of the elastic strain tensor, indicating that {eq:2.5} is isotropic. The symbol $\delta$ denotes the three-dimensional identity tensor. We note that the functional forms of $\lambda$ and $\mu$ in {eq:2.5} can not be arbitrarily chosen; they must derive from an energy functional that is convex in $e$. The damage parameter $\alpha$ is assumed to be governed by an evolution equation of the reaction-diffusion type

$$\tag{eq:2.6}
	\partial_t \alpha = \nabla \cdot \left[ D_\mathrm{a}(\alpha)\nabla \alpha \right] + R_\mathrm{a}(\jmath e, \alpha),
$$

with nonlinear diffusivity $D_\mathrm{a}(\alpha)$ and the reaction/coupling terms collected in $R_\mathrm{a}$.

We revisit the Cauchy momentum balance equation {eq:2.1}, which we combine with {eq:2.5} to yield the concise representation

$$\tag{eq:2.7}
	\partial_t v = r^{-1}\nabla \cdot s(e,\alpha) + R_{\mathrm{v}}(\nabla \alpha, \nabla \nabla \alpha),
$$

with coupling term $R_{v}$. [TODO: density?]

The anelastic part $\bar{e}$ of the strain tensor $\varepsilon$ is governed by Koiter's rule for non-associated plasticity [Koiter, 1953]:

$$\tag{eq:2.8}
	\partial_t \bar{e} = \gamma(\beta) \partial_s G(\jmath s),
$$

in which the function $G$ takes the role of a so-called _plastic potential_, and the scalar coefficient $\gamma = \gamma(\beta)$ is called the _plastic multiplier_. In traditional plasticity models, $\gamma$ is used as a Lagrange multiplier that optimizes [TODO: finish], but here we assume that it is a direct function of a further damage field $\beta = \beta(t,x)$, which itself smoothly evolves according to the reaction-diffusion equation

$$\tag{eq:2.9}
	\partial_t \beta = \nabla \cdot \left[ D_\mathrm{b}(\beta) \nabla \beta \right] + R_\mathrm{b}(\jmath e, \alpha, \beta),
$$

with nonlinear diffusivity $D_\mathrm{b} = D_\mathrm{b}(\beta)$ and the reaction (or source) terms collected in $R_\mathrm{b}$.

Koiter's rule {eq:2.8} can also be directly stated in terms of unknowns $\alpha$, $\beta$, and the invariants of $e$. We make use of this fact, the relation {eq:2.4}, and the differential relation $v = \partial_t u$, to write an evolution equation for the elastic part of the strain rate,

$$\tag{eq:2.10}
	\partial_t e = \nabla^\mathrm{s} v - \partial_t \bar{e} = \nabla^\mathrm{s} v - \gamma(\beta) \partial_s G(\jmath s(e,\alpha)),
$$

Special cases of the rheology {eq:2.4}--{eq:2.7} include the damage-breakage rheology of [Lyakhovski, Ben-Zion, et al., 2011, 2014a, 2014b, etc.], the rate and state rheology of [Pranger et al., 2022], and the temperature and grain-size dependent rheology of [The Geodynamicists]. We will use [one or the other or both] to test our algorithm in [section X]. Ultimately, we seek to find a sufficiently accurate numerical approximation to {eq:2.1}--{eq:2.7} in terms of trajectories of the state fields $u$, $v$, $\bar{e}$, $\alpha$, and $\beta$ in the interior of the domain $\mathrm{T} \otimes \Omega$, given sufficient data on $\partial \mathrm{T}$ and $\partial \Omega$.

> TODO: include temperature

3: Numerical methods
-----------------

Summarizing the previous section, we consider the numerical solution of the system

$$\tag{eq:3.1a}
	\partial_t v = r^{-1}\nabla \cdot s(e,\alpha) + R_{\mathrm{v}}(\nabla \alpha, \nabla \nabla \alpha),
$$

$$\tag{eq:3.1b}
	\partial_t e = \nabla^\mathrm{s} v - \gamma(\beta) \partial_s G(s(\alpha, \jmath e)),
$$

$$\tag{eq:3.1c}
	\partial_t \alpha = \nabla \cdot \left[ D_\mathrm{a}(\alpha)\nabla \alpha \right] + R_\mathrm{a}(\jmath e, \alpha),
$$

$$\tag{eq:3.1d}
	\partial_t \beta = \nabla \cdot \left[ D_\mathrm{b}(\beta) \nabla \beta \right] + R_\mathrm{b}(\jmath e, \alpha, \beta),
$$

which we summarize as the ordinary differential equation (ODE)

$$\tag{eq:3.2}
	\frac{\partial U}{\partial t} = F(U; t),
$$

with $U$ a the vector of unknowns $U = [v, e, \alpha, \beta]^T$ and $F$ a (partial differential) operator that acts on $U$ at a time $t$. Explicit time-dependence is included in case time-varying boundary conditions are needed.

**3.1: Considerations for time discretization**

We use the method of lines to separate time and space discretization of these equations. For the sake of simplicity of implementation we favor second-order accurate methods in both cases, which has the advantage that the combined accuracy of space and time discretization is also close to second-order [citation needed].

When discretizing a system of ODEs like {eq:3.2} in time, we typically have a choice between so-called explicit and implicit methods. Explicit methods are algorithmically simple and have a low, constant computational cost per time step, but the step size is limited to some fraction of the characteristic time scale of the fastest-growing perturbations admitted by the ODE. Typically, perturbations grow much faster than the actual solution of the initial value problem, a phenomenon called stiffness [citation needed]. This problem is further exacerbated by PDE components in the right-hand-side of {eq:3.2}. After spatial discretization [Section 3.6], these admit perturbations whose characteristic time scale grows with the grid distance; linearly for hyperbolic PDE components and even quadratically for parabolic PDE components [citation needed].

Implicit methods on the other hand, do not suffer from time step restrictions when they are A-stable [Section 3.3], but require the solution of a (large) system of simultaneous algebraic equations, which incurs a considerable, non-constant cost per time step [citation needed]. Given a solution method whose algorithmic complexity scales sub-linearly with time step size, implicit schemes can easily outperform explicit schemes whenever the solution is 'slow'[^1] [citation needed]. However, as will become clear by [Section 3.6], our favored class of iterative, matrix-free solvers places various constraints on the properties of the operator $f$ summarized in {eq:3.xyzz}, such as positive (negitive) definiteness or symmetry. The physical processes modeled in {eq:3.1a}--{eq:3.1d} unfortunately do not result in such properties, except in limit cases.

One of the major aims of this paper is the resolution of this problem. We propose an additive decomposition of the governing equations {eq:3.1a}--{eq:3.1d}, in which some of the PDE terms can be decoupled and solved with an implicit method in a routine way, and the remainder with an explicit method using an additive Runge-Kutta (ARK) scheme. We derive subsequently the properties of this time integration scheme, and show that the resulting time step constraints are independent of the spatial resolution, despite the fact that not all of the PDE components are treated implicitly. Finally, we embed into this scheme an error estimator and propose time step controls to keep the discretization error bounded.

[^1]: Here 'slow' is relative: the solution to a friction problem may e.g. be 'slow' both during the 'stick' and the 'slip' phases, and only be 'fast' during the transitions between the two.


**3.2: An additive Runge-Kutta method**

We first consider the additive decomposition and time integration of the generic system {eq:3.2}, and subsequently propose in [Section 3.4] a suitable decomposition of the governing equations {eq:3.1a}--{eq:3.1d}. In this section we follow the work of [Giraldo et al., 2013], who applied the method to atmospheric flow. We symbolically write the additive decomposition of {eq:3.2} into explicitly and implicitly solved terms as

$$\tag{eq:3.3}
	\partial_t U = F(U; t) = F^\mathrm{im}(U; t) + F^\mathrm{ex}(U; t),
$$

in which $F^\mathrm{im}(U; t)$ contains the terms to be integrated implicitly, and $F^\mathrm{ex}(U; t)$ the terms to be integrated explicitly.

We adopt the additive Runge-Kutta scheme

<!-- $$\tag{eq:3.tta}
	w_i  =    u^n + \tau \sum\limits_{j=1}^{i-1} a_{ij} g(w_j) \\
	\hspace{11em} + \tau \sum\limits_{j=1}^{i} \tilde{a}_{ij} \tilde{g}(w_j), \qquad i = 1, \ldots, 3
$$ -->

$$\tag{eq:3.4a}
	\bar{T} = t^n \bar{\mathrm{1}} + \tau\, \bar{C},
$$

$$\tag{eq:3.4b}
	\bar{W} = U^n \bar{\mathrm{1}} + \tau \left[\bar{\mathbf{A}}^\mathrm{im} \bar{F}^\mathrm{im}(\bar{W}; \bar{T}) + \bar{\mathbf{A}}^\mathrm{ex} \bar{F}^\mathrm{ex}(\bar{W}; \bar{T})\right],
$$

$$\tag{eq:3.4c}
	U^{n+1} = U^n + \tau\, \bar{B}^\mathrm{T} \left[\bar{F}^\mathrm{im}(\bar{W}; \bar{T}) + \bar{F}^\mathrm{ex}(\bar{W}; \bar{T}) \right],
$$

in which $n$ denotes the step index, $\tau$ the step size, $t^n = t$, $t^{n+1} = t + \tau$, and $U^n = U(t)$, $U^{n+1} = U(t + \tau)$. The barred symbols denote vectors and square matrices of size $m$, the number of stages in the Runge-Kutta scheme. We use $m = 3$ stages to generate a method that is of second-order accuracy. The intermediate stage values of the solution $U$ are stored in $\bar{W}$, the intermediate stage values of time are stored in $\bar{T}$. The symbol $\bar{1}$ denotes a vector of ones, so that e.g. $t^n \bar{1} = [t^n, t^n, t^n]^\mathrm{T}$. $\bar{\mathbf{A}}^\mathrm{im,ex}$, $\bar{B}$ and $\bar{C}$ denotes coefficient matrices and vectors, subject to the condition that

$$\tag{eq:3.5a}
	\bar{C} = \bar{\mathbf{A}}^\mathrm{im}\bar{1} = \bar{\mathbf{A}}^\mathrm{ex} \bar{1}.
$$

$$\tag{eq:3.5b}
	\bar{B} \cdot \bar{1} = 1.
$$

Finally, the notation $\bar{F}^\mathrm{im,ex}(\bar{W}; \bar{T})$ must be understood to mean

$$
	\bar{F}^\mathrm{im,ex}(\bar{W}; \bar{T}) = \left[F^\mathrm{im,ex}(W_1; T_1),\; F^\mathrm{im,ex}(W_2; T_2),\; F^\mathrm{im,ex}(W_3; T_3)\right]^\mathrm{T}.
$$

Equation {eq:3.4b} defines a system of equations in $m$ unknowns, which can be solved sequentially if the coefficient matrices $\bar{\mathbf{A}}^\mathrm{im,ex}$ are lower diagonal. In this case the scheme is said to be _diagonally implicit_. The function $F^\mathrm{ex}$ is only treated explicitly when the coefficient matrix $\bar{\mathbf{A}}^\mathrm{ex}$ is _strictly_ lower diagonal.

We populate $\bar{\mathbf{A}}^\mathrm{im}$ with coefficients from the TR-BDF2 time integration scheme [CITE, Giraldo 2013], which can be written as an explicit first stage $W_1 = U^{n}$ at $T_1 = t^n$, a second stage at $T_2 = t^n + c_2 \tau$ that is generated by the trapezoidal rule $W_2 = U^{n} + (c_2/2) F^\mathrm{im}(W_1; T_1) + (c_2/2) F^\mathrm{im}(W_2; T_2))$, and a final stage at $T_3 = t^n + \tau$ given by the second-order backwards difference formula (BDF2). The TR-BDF2 scheme would be finished at $U^{n+1} = W_3$, which can be extended into the Runge-Kutta framework by letting $\bar{B}$ be equal to the last row of $\bar{\mathbf{A}}^\mathrm{im}$, i.e.

$$\tag{eq:3.6}
	\bar{B} = \bar{\mathbf{A}}^\mathrm{im} [0, 0, 1]^\mathrm{T}.
$$

However, by keeping the values of $\bar{B}$ open for the moment (subject to {eq:3.5b}), we will be able to generate schemes with different properties at the expense of only one explicit finishing stage compared to the TR-BDF2 scheme (see Section 3.3).

[TODO: Laudable properties of TR-BDF2].

The resulting in the Butcher tableau [CITE] is given by

$$\tag{eq:3.7a}
    \begin{array}{c|ccc}
			  \bar{C}
			& \bar{\mathbf{A}}^\mathrm{im} \\[.7em] \hline\\[-.5em]
			& \bar{B}^\mathrm{T}
		\end{array} \quad = \quad \begin{array}{c|ccc}
			  0
			&  
			&  
			&   \\[.7em]
			  c_2
			& \frac{1}{2}c_2
			& \frac{1}{2}c_2
			&   \\[.7em]
			  1
			& \frac{1}{2}(2-c_2)^{-1}
			& \frac{1}{2}(2-c_2)^{-1}
			& (1-c_2)(2-c_2)^{-1} \\[.7em]\hline\\[-.5em]
			& b_1
			& b_2
			& 1 - b_1 - b_2
		\end{array}.
$$

Constrained by {eq:3.5a} and the condition that $\bar{\mathbf{A}}^\mathrm{ex}$ must be strictly lower triangular, the Butcher tableau for the explicit Runge-Kutta scheme reads

$$\tag{eq:3.7b}
    \begin{array}{c|ccc}
			  \bar{C}
			& \bar{\mathbf{A}}^\mathrm{ex} \\[.7em] \hline\\[-.5em]
			& \bar{B}^\mathrm{T}
		\end{array} \quad = \quad \begin{array}{c|ccc}
			  0
			&  
			&  
			&   \\[.7em]
			  c_2
			& c_2
			& 
			&   \\[.7em]
			  1
			& 1-a^{ex}_{32}
			& a^{ex}_{32}
			& \\[.7em]\hline\\[-.5em]
			& b_1
			& b_2
			& 1 - b_1 - b_2
		\end{array}.
$$

We thus have four parameters, $b_1$, $b_2$, $c_2$ and $a^{ex}_{32}$, that can be tuned to some advantage.

**3.3 Stability and Accuracy**

We analyze the numerical stability and accuracy of the scheme {eq:3.4a--c}, {eq:3.7a}, {eq:3.7b} using Dahlquist's problem [Dahlquist, 1963], which we adapt to our additive IMEX problem as follows:

$$\tag{eq:3.8}
	\tau\, \partial_t U = \zeta U = \zeta^\mathrm{im} U + \zeta^\mathrm{ex} U.
$$

Here, $\zeta^\mathrm{im}/\tau$ and $\zeta^\mathrm{ex}/\tau$ can be interpreted as eigenvalues of the linearization of $F^\mathrm{im}(U; t)$ and $F^\mathrm{ex}(U; t)$ with respect to $U$. We substitute {eq:3.7} into {eq:3.4a--c} to obtain

$$\tag{eq:3.9a}
	\bar{W} = U^n \bar{\mathrm{1}} + \left[\zeta^\mathrm{im} \bar{\mathbf{A}}^\mathrm{im} + \zeta^\mathrm{ex}\bar{\mathbf{A}}^\mathrm{ex} \right] \bar{W},
$$

$$\tag{eq:3.9b}
	U^{n+1} = U^n + \zeta\, \bar{B}^\mathrm{T} \bar{W},
$$

Next, we invert {eq:3.9a} for $\bar{W}$ and substitute the result into {eq:3.9b};

$$\tag{eq:3.10a}
	U^{n+1} = P(\zeta^\mathrm{im}, \zeta^\mathrm{ex}) U^n.
$$

with

$$\tag{eq:3.10b}
	P(\zeta^\mathrm{im}, \zeta^\mathrm{ex}) = 1 + \zeta\, \bar{B}^\mathrm{T} \left[ \bar{\mathbf{I}} - \zeta^\mathrm{im} \bar{\mathbf{A}}^\mathrm{im} - \zeta^\mathrm{ex} \bar{\mathbf{A}}^\mathrm{ex} \right]^{-1} \bar{\mathrm{1}}
$$

L-stability of the implicit end-member scheme ($\zeta^\mathrm{ex} = 0$) is demonstrated by the fact that $\lim\limits_{|\zeta^\mathrm{im}| \to \infty}P(\zeta^\mathrm{im}, 0) = 0$ [e.g. Hosea & Shampine, 1996]. This L-stable property means that the method arrives at $U^{n+1} = 0$ in one step in the limit of infinite time step. For $\zeta^\mathrm{im} \in \mathbb{R}^-$ this equates to the analytical solution 

$$\tag{eq:3.11}
	\tilde{U}^{n+1} = U^{n} \exp{\zeta}
$$

to test problem {eq:3.8}.

<!-- By using $\mathrm{det}\left(\mathbf{M}\right) \mathbf{M}^{-1} = \mathrm{adj}\left(\mathbf{M}\right)$ for any invertible matrix $\mathbf{M}$, {eq:3.10} is concisely stated as

$$\tag{eq:3.11a}
	P^\mathrm{im}(\zeta^\mathrm{im}) U^{n+1} = P^\mathrm{ex}(\zeta^\mathrm{im}, \zeta^\mathrm{ex}) U^n,
$$

in which the (bivariate) polynomials $P^\mathrm{im}(\zeta^\mathrm{im})$ and $P^\mathrm{ex}(\zeta^\mathrm{im}, \zeta^\mathrm{ex})$ can be expressed as

$$\tag{eq:3.11b}
	P^\mathrm{im}(\zeta^\mathrm{im}) = \mathrm{det}\; \left[\bar{\mathbf{I}} - \zeta^\mathrm{im} \bar{\mathbf{A}}^\mathrm{im} - \zeta^\mathrm{ex} \bar{\mathbf{A}}^\mathrm{ex} \right] = \prod \left[\bar{1} - \zeta^\mathrm{im} \mathrm{diag}( \bar{\mathbf{A}}^\mathrm{im})\right],
$$

$$\tag{eq:3.11c}
	P^\mathrm{ex}(\zeta^\mathrm{im}, \zeta^\mathrm{ex}) = P^\mathrm{im}(\zeta^\mathrm{im}) + \zeta\, \bar{B}^\mathrm{T} \left(\mathrm{adj}\; \left[\bar{\mathbf{I}} - \zeta^\mathrm{im} \bar{\mathbf{A}}^\mathrm{im} - \zeta^\mathrm{ex} \bar{\mathbf{A}}^\mathrm{ex} \right] \right) \bar{\mathrm{1}}.
$$

The second equality in {eq:3.11b} is due to the determinant of any lower-triangular matrix being equal to the product of its diagonal.
$P^\mathrm{ex}(\zeta^\mathrm{im}, \zeta^\mathrm{ex})$ is a bivariate polynomial in $\zeta^\mathrm{im}$ and $\zeta^\mathrm{ex}$ in part because an adjugate is ultimately an entry-wise polynomial in its argument. -->

This can be expanded using a Taylor series around the origin, giving

$$\tag{eq:3.12a}
	\tilde{U}^{n+1} = Q(\zeta) U^{n},
$$

with

$$\tag{eq:3.12b}
	Q(\zeta) = \left( 1 + \zeta + \tfrac{1}{2} \zeta^2 + \tfrac{1}{6} \zeta^3 + \mathcal{O}(\zeta^4) \right).
$$

We can now combine {eq:3.12a,b} and {eq:3.11c} and write for the integration error $\tilde{\varepsilon}^{n+1} = U^{n+1} - \tilde{U}^{n+1}$:


$$\tag{eq:3.13}
	\tilde{\varepsilon}^{n+1} = \left[P(\zeta^\mathrm{im}, \zeta^\mathrm{ex}) - Q(\zeta)\right] U^{n},
$$

The polynomial $R(\zeta^\mathrm{im}, \zeta^\mathrm{ex})$ evaluates to

$$\tag{eq:3.14}
	R(\zeta^\mathrm{im}, \zeta^\mathrm{ex}) =
	  (b_1 + b_2 - b_2 c_2 - \tfrac{1}{2}) \mathcal{O}(\zeta^2),
$$

which means the second-order error can be eliminated by choosing

$$\tag{eq:3.15}
	b_1 = \tfrac{1}{2} + b_2 (c_2 - 1).
$$

At this point the second-order accuracy of the TR-BDF2 scheme can be verified, since condition {eq:3.15} is consistent with {eq:3.6}. The remainder after elimination of $b_1$ using {eq:3.15} is

$$\tag{eq:3.16}
	R(\zeta^\mathrm{im}, \zeta^\mathrm{ex}) = (1 + 6 b_2 c_2 (c_2 - 1)) \mathcal{O}(\zeta^3_\mathrm{im}) \\[.8em]
	  + (2(a_{32}^\mathrm{ex}-1) b_2 c_2 - a_{32}^\mathrm{ex} + 2 b_2) \mathcal{O}(\zeta^2_\mathrm{im}\zeta_\mathrm{ex}) \\[.8em]
	  + (2 c_2 ((4 a_{32}^\mathrm{ex}-1) b_2 c_2-2 a_{32}+b_2)+1) \mathcal{O}(\zeta_\mathrm{im}\zeta^2_\mathrm{ex}) \\[.8em]
	  + (3 a_{32}^\mathrm{ex} c_2 \left(2 b_2 c_2-1\right)+1) \mathcal{O}(\zeta^3_\mathrm{ex})
$$

These leading error terms can be eliminated by choosing

$$\tag{eq:3.17}
	a_{32}^\mathrm{ex} = (c_2 - 1)(3 c_2^2 - 2 c_2)^{-1}, \qquad  b_2 = \tfrac{1}{6}(c_2-c_2^2)^{-1}
$$

Finally, the coefficients of the remaining error terms are minimized by choosing $c_2 = 2 - \sqrt{2}$, which coincidentally yields a _singly diagonally implicit_ Runge-Kutta (SDIRK) scheme [CITE], i.e. $a^\mathrm{im}_{22} = a^\mathrm{im}_{33}$. For linear problems solved using assembled Jacobians this property would be advantageous, but this does not apply here. Summarizing, we have the schemes

$$
	\bar{\mathbf{A}}^\mathrm{im} = \left(\begin{array}{ccc}
		0 & 0 & 0 \\[.7em]
		1-\frac{1}{2}\sqrt{2} & 1-\frac{1}{2}\sqrt{2} & 0 \\[.7em]
		\frac{1}{4}\sqrt{2} & \frac{1}{4}\sqrt{2} & 1-\frac{1}{2}\sqrt{2}
	\end{array}\right)
$$

$$
	\bar{\mathbf{A}}^\mathrm{ex} = \left(\begin{array}{ccc}
		0 & 0 & 0 \\[.7em]
		2-\sqrt{2} & 0 & 0 \\[.7em]
		-\frac{1}{2}-\sqrt{2} & \frac{3}{2}+\sqrt{2} & 0
	\end{array}\right)
$$

$$
	\bar{B}_\mathrm{2} = \left(\begin{array}{c}
		\frac{1}{4}\sqrt{2} \\[.7em]
		\frac{1}{4}\sqrt{2} \\[.7em]
		1-\frac{1}{2}\sqrt{2}
	\end{array}\right)
$$

$$
	\bar{B}_\mathrm{3} = \left(\begin{array}{c}
		\frac{1}{12} \left(4 -   \sqrt{2} \right) \\[.7em]
		\frac{1}{12} \left(4 + 3 \sqrt{2} \right) \\[.7em]
		\frac{1}{6 } \left(2 -   \sqrt{2} \right)
	\end{array}\right)
$$



<!--
- three-stage Butcher-Chen method [Butcher and Chen, 2000]
- In line with [Giraldo et al., 2013], w
- For $\gamma = 2-\sqrt{2} \approx 0.58...$, a singly diagonally implicit Runge-Kutta (SDIRK) scheme is obtained for the implicit integrator, meaning that any assembled implicit operator can be reused for both stages. We use matrix-free solvers (Section [TODO]), so we merely gain an aesthetic benefit.
- When the coefficients $b^\mathrm{T} = [b_1, b_2, b_3]$, $\lVert b \rVert_1 = b_1 + b_2 + b_3 = 1$ in {eq:3.13d} are chosen equal to the corresponding coefficients $\tilde{a}_3^\mathrm{T} = [\tilde{a}_{31}, \tilde{a}_{32}, \tilde{a}_{33}]$ in {eq:3.13c}, we obtain a redundant reformulation of the original TR-BDF2 scheme with. The more general form {eq:3.13} allows us however to choose different values of $b^\mathrm{T}$ that also eliminate the third-order truncation error $\mathcal{O}(h_t^3)$, at the expense of L and A stability. Inserting again the test equation $\partial y/\partial t = \lambda y$ and following the procedure that led to {eq:3.11}, we now obtain
- We can see that the scheme is L-stable as long as $\gamma \neq 1$, since then the polynomial order of $P_\mathrm{im} > P_\mathrm{ex}$ and thus $\lim\limits_{h_t \lambda \to\pm\infty} P_\mathrm{ex}(h_t \lambda)/P_\mathrm{im}(h_t \lambda) = 0$ (i.e. for negative $\lambda$, the solution $y$ has a steady state on which it converges as the step size increases).
- The free parameter $\gamma$ is chosen such that for the system of ODEs $\partial y / \partial t = J y$, the operators $P_\mathrm{im,2}(h_t J)$ and $P_\mathrm{im,1}(h_t J)$ are the same and are assembled only once per time step[^1]. Thus,
$$\tag{eq:3.12}
	\frac{\gamma}{2} = \frac{1-\gamma}{2-\gamma} \implies \gamma = 2 \pm \sqrt{2}.
$$
The solution with minus sign, $2 - \sqrt{2} \approx 0.586$ not only lies in the desired interval $(0,1]$, but also yields a much lower coefficient of the $\mathcal{O}(h_t^3)$ truncation error ($\sim 0.04$ vs. $\sim 1.4$) and is therefore selected.
- We implement matrix-free algorithms and so the assembly cost is no issue for us.
-->

**3.4: Error estimation and time step control**

[IGNORE BELOW, TO BE REORGANIZED!]
The more general form {eq:3.13} allows us however to choose different values of $b^\mathrm{T}$ that also eliminate the third-order truncation error $\mathcal{O}(h_t^3)$, at the expense of L and A stability. Inserting again the test equation $\partial y/\partial t = \lambda y$ and following the procedure that led to {eq:3.11}, we now obtain


$$\tag{eq:3.15}
	P^\ast_\mathrm{im}(h_t \lambda) \left(1 + h_t \lambda + \frac{1}{2}(h_t \lambda)^2 + \frac{1}{6}(h_t \lambda)^3 + \mathcal{O}(h_t^4) \right) - P^\ast_\mathrm{ex}(h_t \lambda) = \\[1em]
    = g_1(b) (h_t \lambda) + g_2(b) (h_t \lambda)^2 + g_3(b) (h_t \lambda)^3 + \mathcal{O}(h_t^4),
$$

with the appropriate polynomials $P^\ast_\mathrm{im}$ and $P^\ast_\mathrm{ex}$, and coefficients $g(b)^\mathrm{T} = [g_1(b), g_2(b), g_3(b)]$, which are equated to zero to yield:

$$\tag{eq:3.16}
	P^\ast_\mathrm{im}(h_t \lambda) \left(1 + h_t \lambda + \frac{1}{2}(h_t \lambda)^2 + \frac{1}{6}(h_t \lambda)^3 + \mathcal{O}(h_t^4) \right) - P^\ast_\mathrm{ex}(h_t \lambda) =  \mathcal{O}(h_t^4), \\[1em]
	b_1 = \tfrac{1}{3} - \tfrac{1}{12}\sqrt{2}, \\[1em]
	b_2 = \tfrac{1}{3} + \tfrac{1}{4} \sqrt{2}, \\[1em]
	b_3 = \tfrac{1}{3} - \tfrac{1}{6} \sqrt{2}.
$$

Having computed a stable second-order accurate approximation $y_{n+1}$ to the solution $y(t_n + h_t)$, and a third-order accurate approximation $y^\ast_{n+1}$ to the same, we may estimate the $\mathcal{O}(h^3)$ truncation error $\epsilon_{n+1}$ as
$$\tag{eq:3.17}
	\epsilon_{n+1} - \left(1-\frac{1}{\sqrt{2}}\right) h_t f(\epsilon_{n+1}) = y^\ast_{n+1} - y_{n+1}.
$$
The implicit filter is due to [Shampine 1984, Hosea \& Shampine, 1996; in Bonaventura 2018] and gives the error measure $\epsilon^{n+1}$ the correct limit behavior for large $h_t$.

Finally, we seek a three-stage explicit Runge-Kutta scheme that complements the TR-BDF2 scheme by utilizing the same fractional step $\gamma h_t = (2 - \sqrt{2}) h_t$ and the same final stage coefficients $b^\mathrm{T} = [1/(2\sqrt{2}), 1/(2\sqrt{2}), 1 - 1/\sqrt{2}]$ (compare {eq:3.13}, {eq:3.18}; see [Giraldo, 2013]). We thus seek the coefficient $a_{32}$ in

$$\tag{eq:3.18}
	w_1 = y_n \\[1em]
	w_2 = y_n + h_t \left(2-\sqrt{2}\right) f(w_1), \\[1em]
	w_3 = y_n + h_t \left[ (1 - a_{32}) f(w_1) + a_{32} f(w_2) \right], \\[1em]
	y^\ast_{n+1} = y_n + h_t \left[ \left(\frac{1}{2\sqrt{2}}\right) f(w_1) + \left(\frac{1}{2\sqrt{2}}\right) f(w_2) + \left(1 - \frac{1}{\sqrt{2}}\right) f(w_3) \right],
$$

that maximizes the accuracy of the method. We again construct the update for the test equation $\partial y/\partial t = \lambda y$, $\lambda \in \mathbb{C}$,

$$\tag{eq:3.19a} %\label{eq:expoly}
	y_{n+1} = P(h_t \lambda) y_n,
$$
$$\tag{eq:3.19b}
	P(z) = 1 + z + \frac{1}{2}z^2 + (3 - 2\sqrt{2}) a_{32} z^3
$$
$$\tag{eq:3.19c}
	P(h_t \lambda) - \exp(h_t \lambda) = \left((3 - 2\sqrt{2}) a_{32} - \frac{1}{6} \right) (h_t \lambda)^3 + \mathcal{O}(h_t^4)
$$
$$\tag{eq:3.19d}
	= \mathcal{O}(h_t^4) \iff a_{32} = \tfrac{1}{2} + \tfrac{1}{3}\sqrt{2}.
$$

----

![fig:stab](https://media.githubusercontent.com/media/cpranger/DamageBreakage/main/paper/figures/stability.png)

> Stability regions $\lvert P(h_t \lambda) \rvert \leq 1$ in the complex plane $h_t \lambda \in \mathbb{C}$ derived from the scalar test equation $\partial y / \partial t = \lambda y$, of (top left) the second-order L-stable TR-BDF2 scheme, (top right) its third-order diagonally implicit Runge-Kutta extension, (bottom left) an embedded second-order explicit Runge-Kutta scheme, and (bottom right) an embedded third-order explicit Runge-Kutta scheme. In the latter, a circular region $D$ is plotted with origin $0$ and radius $\sqrt{3}$, the intersection of which with the left half plane is a proper subset of the stability region $\lvert P(h_t \lambda) \rvert \leq 1$ of the explicit 3rd-order Runge-Kutta scheme by which the explicit part of the equations is integrated.

The stability regions of the four schemes derived here (second-order L-stable TR-BDF2, its third-order diagonally implicit Runge-Kutta extension, its associated third-order explicit Runge-Kutta counterpart, and its embedded second-order explicit RK method) are plotted in Figure {fig:stab}. In the bottom-right panel of the same figure, the region of stability of the explicit RK scheme is approximated with a circular region with radius $\sqrt{3}$ centered around the origin. The intersection of this region with the left half complex plane is a proper subset of the complete stability region, and allows us to determine an appropriate time step using only the the spectral radius of the explicit Jacobian $\mathbf{J}_\mathrm{ex}$ (see Section 3.1):

$$\tag{eq:3.20}
	h_t = \sqrt{3}/\rho(\mathbf{J}_\mathrm{ex}).
$$

Since the circular region intersects the imaginary axis at critical stability and passes close by the three complex zeroes of the stability polynomial, this choice gives both excellent conservation properties of the purely hyperbolic components, and excellent damping properties of the stiff dissipative components. We will refine the issue of time step selection in Section 3.4.


Inspired by the multi-rate extension of TR-BDF2 of [Bonaventura et al., 2018], we design the following much simplified mechanism for controlling the error on the implicitly solved components of {eq:3.4}. We estimate the implicit error using {eq:3.21c}. Then, given the absolute and relative tolerances $\tau_\mathrm{a}$ and $\tau_\mathrm{r}$, we compute the dimensionless error measure

$$\tag{eq:3.23}
	\eta_i = \frac{\epsilon_{n+1,i}}{\tau_\mathrm{r} \lvert y_{n+1,i} \rvert + \tau_\mathrm{a}}.
$$

We then compute the $L_\infty$ norm of $\eta$ over the vector blocks $[\vec{v}_{n+1},\mathbf{e}_{n+1}]$ (since these two will be solved together in block reduced form), $[\alpha_{n+1}]$, and $[\beta_{n+1}]$. We denote these maximum dimensionless errors by $\eta_\mathrm{v}$, $\eta_\mathrm{a}$, and $\eta_\mathrm{b}$, respectively.

If any of these measures exceed a value of one, a new time step is computed by [Bonaventura 2018 and references [14], [19] therein]

$$\tag{eq:3.24}
	h_t^\ast = \nu h_t \mathrm{max}(\eta_\mathrm{v}, \eta_\mathrm{a}, \eta_\mathrm{b})^{-1/3},
$$

with $\nu \in (0, 1)$ a safety coefficient and the power $1/3$ specific to a method with second-order accuracy.

The time step is then redone with $h_t^\ast$ and this process is repeated until $\mathrm{max}(\eta_\mathrm{v}, \eta_\mathrm{a}, \eta_\mathrm{b}) \leq 1$. When this condition is satisfied, there ought to be some lingering memory of the failure of the explicit step, and the time step selection {eq:3.20} is best amended with the recursion

$$\tag{eq:3.25} %\label{eq:timestep2}
	h_t^{(n+1)} = \mathrm{min}( \sqrt{3}/\rho(\mathbf{J}_\mathrm{ex}), (1/\nu) h_t^{(n)} ),
$$

with the same safety coefficient $\nu \in (0, 1)$.

At each stage of time step refinement the initial guess of the new solution may be significantly improved by exploiting the globally $C_1$-continuous cubic Hermite interpolation that is given in [Section 5 of Bonaventura et al., 2018], which requires only the known solution stages in its evaluation.

[It might actually be possible to save some refinements by solving a constraint optimization problem using the cubic Hermite interpolation of the solution over the time step, optimizing the non-dimensional error measure $\eta$ subject to the constraint $\eta \geq \nu$.]

<!--
$$\tag{eq:3.21c} %\label{eq:imexarkc}
	(\mathbf{I} - \tilde{b}_{s} h_t \mathbf{J}_\mathrm{im}) \epsilon^{n+1} = h_t \sum\limits_{j=1}^{s} (b^\ast_{j} - b_{j}) \tilde{F}(w_j),
$$


$$\tag{eq:3.22d}
	b^\ast = \left[
		\tfrac{1}{3} - \tfrac{1}{12}\sqrt{2},\;\; \tfrac{1}{3} + \tfrac{1}{4} \sqrt{2},\;\; \tfrac{1}{3} - \tfrac{1}{6} \sqrt{2}
	\right]^\mathrm{T}
$$
 -->


**3.5: Decomposition of the governing equations**

We also consider an approximation to this system:

$$\tag{eq:3.2a}
	\frac{\partial\vec{v}}{\partial t} = r^{-1}\nabla \cdot \left[ \mathbf{S}(\jmath e,\alpha_0) e \right] =  \tilde{F}_{\mathrm{v}}(e),
$$

$$\tag{eq:3.2b}
	\frac{\partial e}{\partial t} = \nabla^\mathrm{s} \vec{v} = \tilde{F}_{\mathrm{e}}(\nabla^\mathrm{s} \vec{v}),
$$

$$\tag{eq:3.2c}
	\frac{\partial\alpha}{\partial t} = \nabla \cdot \left[ D_\alpha(\alpha)\nabla \alpha \right] + R_\alpha( \jmath e_0, \alpha) = \tilde{F}_\mathrm{a}(\alpha),
$$

$$\tag{eq:3.2d}
	\frac{\partial\beta}{\partial t} = \nabla \cdot \left[ D_\beta(\beta) \nabla \beta \right] + R_\beta(\jmath e_0, \alpha_0, \beta) =  \tilde{F}_\mathrm{b}(\beta),
$$

which is evaluated given some known $e_0$, $\alpha_0$, and $\beta_0$. We will ensure that

$$\tag{eq:3.3a}
	\frac{\partial F_{\mathrm{e}}}{\partial \vec{v}} = \frac{\partial \tilde{F}_{\mathrm{e}}}{\partial \vec{v}}
$$

$$\tag{eq:3.3b}
	\frac{\partial F_{\mathrm{a}}}{\partial \nabla \alpha} = \frac{\partial \tilde{F}_{\mathrm{a}}}{\partial \nabla \alpha}
$$

$$\tag{eq:3.3c}
	\frac{\partial F_{\mathrm{b}}}{\partial \nabla \beta} = \frac{\partial \tilde{F}_{\mathrm{b}}}{\partial \nabla \beta}
$$

Using this approximation, the complete system is then again written as

$$\tag{eq:3.4a}
	\frac{\partial\vec{v}}{\partial t} = \left[F_{\mathrm{v}}(e, \alpha) - \tilde{F}_{\mathrm{v}}(e)\right] + \tilde{F}_{\mathrm{v}}(e), 
$$

$$\tag{eq:3.4b}
	\frac{\partial e}{\partial t} = \left[F_{\mathrm{e}}(\nabla^\mathrm{s} \vec{v}, \alpha, \beta, \jmath e) - \tilde{F}_{\mathrm{e}}(\nabla^\mathrm{s} \vec{v})\right] + \tilde{F}_{\mathrm{e}}(\nabla^\mathrm{s} \vec{v}), 
$$

$$\tag{eq:3.4c}
	\frac{\partial\alpha}{\partial t} = \left[F_\mathrm{a}( \jmath e, \alpha ) - \tilde{F}_\mathrm{a}(\alpha)\right] + \tilde{F}_\mathrm{a}(\alpha), 
$$

$$\tag{eq:3.4d}
	\frac{\partial\beta}{\partial t}  = \left[F_\mathrm{b}(\jmath e, \alpha, \beta) - \tilde{F}_\mathrm{b}(\beta)\right] + \tilde{F}_\mathrm{b}(\beta),
$$

We will use an implicit-explicit (IMEX) time integration scheme to solve the terms in square brackets explicit in time, and the remainders implicit in time.

The explicit part of {eq:3.4a}--{eq:3.4d} is assigned a Jacobian $\mathbf{J}_\mathrm{ex}$, and the implicit part a Jacobian $\mathbf{J}_\mathrm{im}$, both of which are given by

$$\tag{eq:3.5a}
	\mathbf{J}_\mathrm{ex} = \begin{bmatrix}
		  0
		& \frac{\partial F_{\mathrm{v}}}{\partial \mathbf{e}} - \frac{\partial \tilde{F}_{\mathrm{v}}}{\partial \mathbf{e}}
		& \frac{\partial F_{\mathrm{v}}}{\partial \alpha}
		& 0 \\[.7em]
		  0 % \left\{ \frac{\partial F_{\mathrm{e}}}{\partial \vec{v}} - \frac{\partial \tilde{F}_{\mathrm{e}}}{\partial \vec{v}} = 0 \right\}
		& \frac{\partial F_{\mathrm{e}}}{\partial \mathbf{e}}
		& \frac{\partial F_{\mathrm{e}}}{\partial \alpha}
		& \frac{\partial F_{\mathrm{e}}}{\partial \beta} \\[.7em]
		  0
		& \frac{\partial F_{\mathrm{a}}}{\partial \mathbf{e}}
		& \frac{\partial F_{\mathrm{a}}}{\partial \alpha} - \frac{\partial \tilde{F}_{\mathrm{a}}}{\partial \alpha}
		& 0 \\[.7em]
		  0
		& \frac{\partial F_{\mathrm{b}}}{\partial \mathbf{e}}
		& \frac{\partial F_{\mathrm{b}}}{\partial \alpha}
		& \frac{\partial F_{\mathrm{b}}}{\partial \beta} - \frac{\partial \tilde{F}_{\mathrm{b}}}{\partial \beta}
	\end{bmatrix},
$$

$$\tag{eq:3.5b}
	\mathbf{J}_\mathrm{im} = \begin{bmatrix}
		  0
		& \frac{\partial \tilde{F}_{\mathrm{v}}}{\partial \mathbf{e}}
		& 0
		& 0 \\[.7em]
		  \frac{\partial \tilde{F}_{\mathrm{e}}}{\partial \vec{v}}
		& 0
		& 0
		& 0 \\[.7em]
		  0
		& 0
		& \frac{\partial \tilde{F}_{\mathrm{a}}}{\partial \alpha}
		& 0 \\[.7em]
		  0
		& 0
		& 0
		& \frac{\partial \tilde{F}_{\mathrm{b}}}{\partial \beta}
	\end{bmatrix}.
$$

<!-- The term in the first column of $\mathbf{J}_\mathrm{ex}$ cancels due to {eq:3.3a}. -->

In order to determine a stable time step [Section X], we would like to determine the spectral radius $\rho(\mathbf{J}_\mathrm{ex})$ of $\mathbf{J}_\mathrm{ex}$, which is the largest element by absolute value of its (complex) spectrum $\sigma(\mathbf{J}_\mathrm{ex})$. We can create the block structure

$$\tag{eq:3.6}
	\mathbf{J}_\mathrm{ex} = \left[\begin{array}{c|ccc}
			  0
			& \frac{\partial F_{\mathrm{v}}}{\partial \mathbf{e}} - \frac{\partial \tilde{F}_{\mathrm{v}}}{\partial \mathbf{e}}
			& \frac{\partial F_{\mathrm{v}}}{\partial \alpha}
			& 0 \\[.7em]\hline\\[-.5em]
			  0
			& \frac{\partial F_{\mathrm{e}}}{\partial \mathbf{e}}
			& \frac{\partial F_{\mathrm{e}}}{\partial \alpha}
			& \frac{\partial F_{\mathrm{e}}}{\partial \beta} \\[.7em]
			  0
			& \frac{\partial F_{\mathrm{a}}}{\partial \mathbf{e}}
			& \frac{\partial F_{\mathrm{a}}}{\partial \alpha} - \frac{\partial \tilde{F}_{\mathrm{a}}}{\partial \alpha}
			& 0 \\[.7em]
			  0
			& \frac{\partial F_{\mathrm{b}}}{\partial \mathbf{e}}
			& \frac{\partial F_{\mathrm{b}}}{\partial \alpha}
			& \frac{\partial F_{\mathrm{b}}}{\partial \beta} - \frac{\partial \tilde{F}_{\mathrm{b}}}{\partial \beta}
		\end{array}\right] = \begin{bmatrix} 0 & \vec{\mathbf{B}}^\mathrm{T} \\[.7em] \vec{\mathbf{0}} & \mathbf{D}\;\, \end{bmatrix},
$$

and use the Schur determinant theorem,

$$\tag{eq:3.7}
	\mathrm{det}\begin{pmatrix} A & B \\ C & D \end{pmatrix} = \mathrm{det}(A) \mathrm{det}(D - C A^{-1} B) = \mathrm{det}(D)\mathrm{det}(A - B D^{-1} C),
$$

so that

$$\\
    \sigma(\mathbf{J}_\mathrm{ex}) = \{ \lambda \in \mathbb{C} : \mathrm{det}( \mathbf{J}_\mathrm{ex} - \lambda \mathbf{I} ) = 0 \}
$$

$$\\
	= \{ \lambda \in \mathbb{C} : \mathrm{det}(-\lambda I_\mathrm{v})\mathrm{det}( \mathbf{D} - \lambda \mathbf{I}_\mathbf{D} - \vec{\mathbf{0}}(-\lambda I_\mathrm{v})^{-1}\vec{\mathbf{B}}^\mathrm{T} ) = 0 \}
$$

$$\\
	= \{ \lambda \in \mathbb{C} : \mathrm{det}(-\lambda I_\mathrm{v})\mathrm{det}( \mathbf{D} - \lambda \mathbf{I}_\mathbf{D} ) = 0 \}
$$

$$\tag{eq:3.8}
	= \sigma(0) \cup \sigma(\mathbf{D}) = \{ 0 \} \cup \sigma(\mathbf{D}).
$$

On account of conditions {eq:3.3b}--{eq:3.3c}, all elements of the matrix $\mathbf{D}$ are diagonal submatrices, and the remaining eigenvalue problem could be solved point-wise and in parallel. The multidiagonal matrices $\frac{\partial F_{\mathrm{v}}}{\partial \mathbf{e}} - \frac{\partial \tilde{F}_{\mathrm{v}}}{\partial \mathbf{e}}$ and $\frac{\partial F_{\mathrm{v}}}{\partial \alpha}$ that are contained in the block $\vec{\mathbf{B}}$ are eliminated from the eigenvalue problem {eq:3.8} due to the Schur determinant theorem {eq:3.7} and the vector of zero blocks $\vec{\mathbf{0}}$.

Should we need the spectrum of the implicitly solved system's Jacobian $\mathbf{J}_\mathrm{im}$ in {eq:3.5b}, we can make use of its partly disjoint block structure in combination with the property

$$\tag{eq:3.u}
    \sigma\begin{pmatrix}
        0 & B \\
        C & 0
    \end{pmatrix} = \pm \sqrt{\sigma\begin{pmatrix}
        0 & B \\
        C & 0
    \end{pmatrix}^2} = \pm \sqrt{\sigma\begin{pmatrix}
        B C & 0 \\
        0 & C B
    \end{pmatrix}}.% = \pm \sqrt{\sigma(B C) \cup \sigma(C B)}
$$

to write

$$\\\tag{eq:3.v}
    \sigma(\mathbf{J}_\mathrm{im}) = \pm \sqrt{\sigma \left(\frac{\partial \tilde{F}_{\mathrm{v}}}{\partial \mathbf{e}} \frac{\partial \tilde{F}_{\mathrm{e}}}{\partial \vec{v}}\right) }
    \cup \pm \sqrt{\sigma \left(\frac{\partial \tilde{F}_{\mathrm{e}}}{\partial \vec{v}} \frac{\partial \tilde{F}_{\mathrm{v}}}{\partial \mathbf{e}}\right) }
    \cup \sigma\left(\frac{\partial \tilde{F}_{\mathrm{a}}}{\partial \alpha}\right)
    \cup \sigma\left(\frac{\partial \tilde{F}_{\mathrm{b}}}{\partial \beta}\right).
$$

The square roots and plusminus signs in {eq:3.u} and {eq:3.v} should be understood to be applied element-wise. Since the operators in {eq:3.v} are all real negative (semi-)definite, the spectrum $\sigma(\mathbf{J}_\mathrm{im})$ is anticipated to occupy an interval on the negative real axis, and an interval on the imaginary axis.
