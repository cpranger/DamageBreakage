Numerical methods for modeling of multi-state rock rheology
=====

Casper Pranger, Dave May, Alice Gabriel, et al.

Abstract
--------
We present a numerical scheme for efficiently modeling the rheological evolution of rocks in response to changes in temperature, damage degree, and grain size. The intended range of application includes evolution towards rheological instability and subsequent recovery in lab- or crustal scale domains of up to three spatial dimensions. To this end, elastodynamic effects are modeled as required by the rheological time scale. We achieve these aims at reasonable computational cost through the use of an additive implicit-explicit (IMEX) split Runge-Kutta scheme with error estimation and step size control, in which individual PDE terms are solved implicitly, and coupling terms are solved explicitly. This results in a method in which numerical stability is not limited by spatial resolution.



1:  Introduction
---------------

2: Problem statement (old text; to be revised)
-----------------

Let $t \in \mathrm{T}$ denote a coordinate in a half-open interval $\mathrm{T} = [0,l_0) \subset \mathbb{R}$ of time with boundary $\partial \mathrm{T} = \{0\}$, and
let $\vec{x} \in \Omega$ be a coordinate in a rectangular region $\Omega = [0,l_1] \otimes [0,l_2] \otimes [0,l_3] \subset \mathbb{R}^3$ of 3-space with boundary $\partial\Omega$.

We wish to examine the behavior of an elastic solid occupying the region $\mathrm{T} \otimes \Omega$ with its initial configuration specified on $\partial \mathrm{T}$, that is steadily loaded by displacement boundary conditions on $\partial\Omega_u \subseteq \partial\Omega$, and occasionally rapidly unloaded by transient anelastic processes, by transient material degradation, or by combinations thereof. To this end we consider Cauchy's linear momentum balance law
$$ %\label{eq:momentum}
	\frac{\tilde{d}\vec{v}}{dt} = \frac{1}{r} \nabla \cdot \boldsymbol{s} + \vec{g},
$$ (2.1)
together with the mass balance law
$$ %\label{eq:mass}
	\frac{\tilde{d}r}{dt} = -r \nabla \cdot \vec{v},
$$ (2.2)
in which $r = r(t,\vec{x})$ denotes the mass density, $\vec{v} = \vec{v}(t,\vec{x})$ the material velocity, $\boldsymbol{s} = \boldsymbol{s}(t,\vec{x})$ Cauchy's stress tensor, and $\vec{g} = \vec{g}(t,\vec{x})$ represents the action of body forces, gravitational or otherwise. Body forces are ignored from here on ($\vec{g}(t,\vec{x}) = 0$) at no loss of generality -- although the physical framework developed here ignores self-gravitation. At considerable loss of generality, but required for the sake of concise argument, is our choice to ignore advection and distortion within the elastic body. The material derivative $\tilde{d}/dt$, that is to be evaluated along the characteristic of flow, is replaced by the partial derivative $\partial/\partial t$. The divergence $\nabla \cdot \boldsymbol{s}$ of the Cauchy stress tensor $\boldsymbol{s}$, which formally is to be evaluated in coordinates $\vec{x} + \vec{u}$ distorted by displacement $\vec{u} = \vec{u}(t,\vec{x})$, is instead evaluated in the reference coordinate system $\vec{x}$. The density deviations due to volume change are assumed to have a negligible impact on the dynamic momentum balance, i.e. $\tilde{d}r/dt = \partial r/\partial t = 0$. The displacement $\vec{u} = \vec{u}(t,\vec{x})$ itself is governed by
$$% \label{eq:displacement}
	\frac{\tilde{d}\vec{u}}{dt} = F_{\vec{u}}(\vec{v}) = \vec{v} \approx \frac{\partial \vec{u}}{\partial t}.
$$ (2.3)
We will revisit the implications of these choices in the context of the present methodological work in the discussion in Section [Discussion].

<!-- %Combining the mass balance law (2.1) with \eqref{eq:displacement} and solving for $r = r(\vec{u})$ gives
% \begin{align} \label{eq:density}
	% r(\vec{u}) = r_0 \exp \nabla \cdot (\vec{u} - \vec{u}_0) = r_0 \exp \nabla \cdot \vec{u},  \quad \vec{u}(t = 0,\vec{x}) &= \vec{u}_0(\vec{x}) = \vec{0}, \\
	% r(t=0,\vec{x}) &= r_0(\vec{x}). \nonumber
% \end{align}
% This solution eliminates differential equation (2.1) and unknown $r$ from the problem. -->

In line with our assumption of a sufficiently distortion-free elastic medium is the use of the (symmetric) infinitesimal strain tensor
$$ %\label{eq:inftalstrain}
	\boldsymbol{\varepsilon} = \nabla^\mathrm{s} \vec{u} = \frac{1}{2} \left[ ( \nabla \vec{u})^\mathrm{T} + (\nabla \vec{u}) \right] = \boldsymbol{e} + \bar{\boldsymbol{e}},
$$ (2.4)
which, as indicated, is additively decomposed into an elastic part $\boldsymbol{e} = \boldsymbol{e}(t,\vec{x})$ and an anelastic part $\bar{\boldsymbol{e}} = \bar{\boldsymbol{e}}(t,\vec{x})$.
<!-- % We choose to invert the relationship (2.4) for the function
% \begin{align} \label{eq:inftalstrain2}
	% \boldsymbol{e}_\mathrm{e} = \boldsymbol{e}_\mathrm{e}(\nabla^\mathrm{s} \vec{u}, \bar{\boldsymbol{e}}) = \nabla^\mathrm{s} \vec{u} - \bar{\boldsymbol{e}}.
% \end{align} -->

An isotropic but otherwise quite general nonlinear elastic constitutive equation with structural stresses [e.g. Lyakhovksy et al., 2014b] is given by
$$ %\label{eq:elasticconst}
	\boldsymbol{s}(\boldsymbol{e}, \alpha) = \lambda(\vec{x},\alpha, \jmath(\boldsymbol{e})) \,\boldsymbol{\delta} \;\mathrm{tr}\; \boldsymbol{e} + 2 \mu(\vec{x},\alpha, \jmath(\boldsymbol{e})) \boldsymbol{e} - \vartheta(\vec{x}) (\nabla \alpha) \otimes (\nabla \alpha).
$$ (2.5)
Here $\lambda$ and $\mu$ are the so-called Lam\'e parameters of an isotropic elastic solid, and $\vartheta$ an additional elastic modulus associated with the structural stresses incurred by the damage field $\alpha = \alpha(t,\vec{x})$. The notation $\jmath(e)$ refers to the three scalar invariants $(J_1, J_2, J_3)$ of the elastic strain tensor, indicating that (2.5) is isotropic. The symbol $\boldsymbol{\delta}$ denotes the three-dimensional identity tensor. We note that the functional forms of $\lambda$ and $\mu$ in (2.5) can not be arbitrarily chosen; they must derive from an energy functional that is convex in $\boldsymbol{e}$. The damage parameter $\alpha$ is assumed to be governed by an evolution equation of the reaction-diffusion type
$$ %\label{eq:alphaode}
	\frac{\partial\alpha}{\partial t} = \nabla \cdot \left[ D_\mathrm{a}(\alpha)\nabla \alpha \right] + R_\mathrm{a}( \jmath(\boldsymbol{e}), \alpha) = F_\mathrm{a}( \jmath(\boldsymbol{e}), \alpha ),
$$ (2.6)
with nonlinear diffusivity $D_\mathrm{a} = D_\mathrm{a}(\alpha)$ and the reaction (or source) terms collected in $R_\mathrm{a}$. Reaction and diffusion are further grouped in the function $F_\mathrm{a}$, which incorporates the differential operator $\nabla$.

We revisit the Cauchy momentum balance equation (2.1), which we combine with (2.5) to yield the concise representation
$$ %\label{eq:momentum2}
	\frac{\partial\vec{v}}{\partial t} = r^{-1}\nabla \cdot \left[ \mathbf{S}(\jmath(\boldsymbol{e}),\alpha) \boldsymbol{e} \right] + R_{\vec{v}}(\alpha) =  F_{\mathrm{v}}(\boldsymbol{e}, \alpha),
$$ (2.7)
with fourth-order stiffness tensor $\mathbf{S}$ and reaction/source term $R_{\vec{v}}$, both of which have access to the differential operator $\nabla$. [REDO] [TODO: density?]

The anelastic part $\bar{\boldsymbol{e}}$ of the strain tensor $\boldsymbol{\varepsilon}$ is governed by Koiter's rule for non-associated plasticity [Koiter, 1953]:
$$
	\frac{\partial \bar{\boldsymbol{e}}}{\partial t} = \gamma(\beta) \frac{\partial }{\partial \boldsymbol{s}} G(\jmath(\boldsymbol{s})),
$$ (2.8)
in which the function $G$ takes the role of a so-called _plastic potential_, and the scalar coefficient $\gamma = \gamma(\beta)$ is called the _plastic multiplier_. In traditional plasticity models, $\gamma$ is used as a Lagrange multiplier that optimizes [TODO: finish], but here we assume that it is a direct function of a further damage field $\beta = \beta(t,\vec{x})$, which itself smoothly evolves according to the reaction-diffusion equation
$$ %\label{eq:betaode}
	\frac{\partial\beta}{\partial t} = \nabla \cdot \left[ D_\beta(\beta) \nabla \beta \right] + R_\beta(\jmath(\boldsymbol{e}), \alpha, \beta) =  F_\mathrm{b}(\jmath(\boldsymbol{e}), \alpha, \beta),
$$ (2.9)
with nonlinear diffusivity $D_\mathrm{b} = D_\mathrm{b}(\beta)$ and the reaction (or source) terms collected in $R_\mathrm{b}$. As in (2.6), reaction and diffusion are further grouped in the function $F_\mathrm{b}$, which incorporates the differential operator $\nabla$.

Koiter's rule (2.8) can also be directly stated in terms of unknowns $\alpha$, $\beta$, and the invariants of $\boldsymbol{e}$. We make use of this fact, and the relation (2.4), to write an evolution equation for the elastic part of the strain rate,
$$ %\label{eq:koiter2}
	\frac{\partial \boldsymbol{e}}{\partial t} = \nabla^\mathrm{s} \vec{v} - \frac{\partial \bar{\boldsymbol{e}}}{\partial t} = \nabla^\mathrm{s} \vec{v} - \gamma(\beta) G'(\alpha, \jmath(\boldsymbol{e})) = F_{\mathrm{e}}(\nabla^\mathrm{s} \vec{v}, \alpha, \beta, \jmath(\boldsymbol{e})),
$$ (2.10)
$$
	G'(\alpha, \jmath(\boldsymbol{e})) = \frac{\partial G(\jmath(\boldsymbol{s}(\boldsymbol{e},\alpha)))}{\partial \boldsymbol{s}(\boldsymbol{e},\alpha)}.
$$ (2.11)

Special cases of the rheology (2.4)--(2.7) include the damage-breakage rheology of [Lyakhovski, Ben-Zion, et al., 2011, 2014a, 2014b, etc.], the rate and state rheology of [Pranger et al., 2022], and the temperature and grain-size dependent rheology of [The Geodynamicists]. We will use [one or the other or both] to test our algorithm in [section X]. Ultimately, we seek to find a sufficiently accurate numerical approximation to (2.1)--(2.7) in terms of trajectories of the state fields $\vec{u}$, $\vec{v}$, $\bar{\boldsymbol{e}}$, $\alpha$, and $\beta$ in the interior of the domain $\mathrm{T} \otimes \Omega$, given sufficient data on $\partial \mathrm{T}$ and $\partial \Omega$.

> TODO: include temperature

3: Numerical methods (old text; to be revised)
-----------------

**3.1: IMEX split**

All in all, we consider the system
$$ %\label{eq:system}
	\frac{\partial\vec{v}}{\partial t} = r^{-1}\nabla \cdot \left[ \mathbf{S}(\jmath(\boldsymbol{e}),\alpha) \boldsymbol{e} \right] + R_{\vec{v}}(\alpha) =  F_{\mathrm{v}}(\boldsymbol{e}, \alpha),
$$ (3.1a)
$$
	\frac{\partial \boldsymbol{e}}{\partial t} = \nabla^\mathrm{s} \vec{v} - \gamma(\beta) G'(\alpha, \jmath(\boldsymbol{e})) = F_{\mathrm{e}}(\nabla^\mathrm{s} \vec{v}, \alpha, \beta, \jmath(\boldsymbol{e})),
$$ (3.1b)
$$
	\frac{\partial\alpha}{\partial t} = \nabla \cdot \left[ D_\alpha(\alpha)\nabla \alpha \right] + R_\alpha( \jmath(\boldsymbol{e}), \alpha) = F_\mathrm{a}( \jmath(\boldsymbol{e}), \alpha ),
$$ (3.1c)
$$
	\frac{\partial\beta}{\partial t} = \nabla \cdot \left[ D_\beta(\beta) \nabla \beta \right] + R_\beta(\jmath(\boldsymbol{e}), \alpha, \beta) =  F_\mathrm{b}(\jmath(\boldsymbol{e}), \alpha, \beta).
$$ (3.1d)
We also consider an approximation to this system:
$$ %\label{eq:approxsystem}
	\frac{\partial\vec{v}}{\partial t} = r^{-1}\nabla \cdot \left[ \mathbf{S}(\jmath(\boldsymbol{e}),\alpha_0) \boldsymbol{e} \right] =  \tilde{F}_{\mathrm{v}}(\boldsymbol{e}),
$$ (3.2a)
$$
	\frac{\partial \boldsymbol{e}}{\partial t} = \nabla^\mathrm{s} \vec{v} = \tilde{F}_{\mathrm{e}}(\nabla^\mathrm{s} \vec{v}),
$$ (3.2b)
$$
	\frac{\partial\alpha}{\partial t} = \nabla \cdot \left[ D_\alpha(\alpha)\nabla \alpha \right] + R_\alpha( \jmath(\boldsymbol{e}_0), \alpha) = \tilde{F}_\mathrm{a}(\alpha),
$$ (3.2c)
$$
	\frac{\partial\beta}{\partial t} = \nabla \cdot \left[ D_\beta(\beta) \nabla \beta \right] + R_\beta(\jmath(\boldsymbol{e}_0), \alpha_0, \beta) =  \tilde{F}_\mathrm{b}(\beta),
$$ (3.2d)
which is evaluated given some known $\boldsymbol{e}_0$, $\alpha_0$, and $\beta_0$. We will assume that
<!-- \label{eq:imexconditions} -->
$$ %\label{eq:imexconditionsa}
	\frac{\partial F_{\mathrm{e}}}{\partial \vec{v}} = \frac{\partial \tilde{F}_{\mathrm{e}}}{\partial \vec{v}} 
$$ (3.3a)
$$ %\label{eq:imexconditionsb}
	\frac{\partial F_{\mathrm{a}}}{\partial \nabla \alpha} = \frac{\partial \tilde{F}_{\mathrm{a}}}{\partial \nabla \alpha} 
$$ (3.3b)
$$ %\label{eq:imexconditionsc}
	\frac{\partial F_{\mathrm{b}}}{\partial \nabla \beta} = \frac{\partial \tilde{F}_{\mathrm{b}}}{\partial \nabla \beta}
$$ (3.3c)
Using this approximation, the complete system is then again written as
$$ %\label{eq:imexsystem}
	\frac{\partial\vec{v}}{\partial t} = \left[F_{\mathrm{v}}(\boldsymbol{e}, \alpha) - \tilde{F}_{\mathrm{v}}(\boldsymbol{e})\right] + \tilde{F}_{\mathrm{v}}(\boldsymbol{e}), 
$$ (3.4a)
$$
	\frac{\partial \boldsymbol{e}}{\partial t} = \left[F_{\mathrm{e}}(\nabla^\mathrm{s} \vec{v}, \alpha, \beta, \jmath(\boldsymbol{e})) - \tilde{F}_{\mathrm{e}}(\nabla^\mathrm{s} \vec{v})\right] + \tilde{F}_{\mathrm{e}}(\nabla^\mathrm{s} \vec{v}), 
$$ (3.4b)
$$
	\frac{\partial\alpha}{\partial t} = \left[F_\mathrm{a}( \jmath(\boldsymbol{e}), \alpha ) - \tilde{F}_\mathrm{a}(\alpha)\right] + \tilde{F}_\mathrm{a}(\alpha), 
$$ (3.4c)
$$
	\frac{\partial\beta}{\partial t}  = \left[F_\mathrm{b}(\jmath(\boldsymbol{e}), \alpha, \beta) - \tilde{F}_\mathrm{b}(\beta)\right] + \tilde{F}_\mathrm{b}(\beta),
$$ (3.4d)
We will use an implicit-explicit (IMEX) time integration scheme to solve the terms in square brackets explicit in time, and the remainders implicit in time.

The explicit part of (3.4a)--(3.4d) is assigned a Jacobian $\mathbf{J}_\mathrm{ex}$, and the implicit part a Jacobian $\mathbf{J}_\mathrm{im}$, both of which are given by
$$
	\mathbf{J}_\mathrm{ex} = \begin{bmatrix}
		  0
		& \frac{\partial F_{\mathrm{v}}}{\partial \mathbf{e}} - \frac{\partial \tilde{F}_{\mathrm{v}}}{\partial \mathbf{e}}
		& \frac{\partial F_{\mathrm{v}}}{\partial \alpha}
		& 0 \\[.7em]
		  \left\{ \frac{\partial F_{\mathrm{e}}}{\partial \vec{v}} - \frac{\partial \tilde{F}_{\mathrm{e}}}{\partial \vec{v}} = 0 \right\}
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
$$ (3.5a)
$$
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
$$ (3.5b)
The term in the first column of $\mathbf{J}_\mathrm{ex}$ cancels due to (3.3a).

In order to determine a stable time step [Section X], we would like to determine the spectral radius $\rho(\mathbf{J}_\mathrm{ex})$ of $\mathbf{J}_\mathrm{ex}$, which is the largest element by absolute value of its (complex) spectrum $\sigma(\mathbf{J}_\mathrm{ex})$. We can create the block structure
$$
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
$$ (3.6)
and use the Schur determinant theorem,
$$ %\label{eq:schurdet}
	\mathrm{det}\begin{pmatrix} A & B \\ C & D \end{pmatrix} = \mathrm{det}(A) \mathrm{det}(D - C A^{-1} B) = \mathrm{det}(D)\mathrm{det}(A - B D^{-1} C),
$$ (3.7)
so that
$$ %\label{eq:jacspectrum}
	\sigma(\mathbf{J}_\mathrm{ex}) = \{ \lambda \in \mathbb{C} : \mathrm{det}( \mathbf{J}_\mathrm{ex} - \lambda \mathbf{I} ) = 0 \}
$$
$$
	= \{ \lambda \in \mathbb{C} : \mathrm{det}(-\lambda I_\mathrm{v})\mathrm{det}( \mathbf{D} - \lambda \mathbf{I}_\mathbf{D} - \vec{\mathbf{0}}(-\lambda I_\mathrm{v})^{-1}\vec{\mathbf{B}}^\mathrm{T} ) = 0 \}
$$
$$
	= \{ \lambda \in \mathbb{C} : \mathrm{det}(-\lambda I_\mathrm{v})\mathrm{det}( \mathbf{D} - \lambda \mathbf{I}_\mathbf{D} ) = 0 \}
$$
$$
	= \sigma(0) \cup \sigma(\mathbf{D}) = \{ 0 \} \cup \sigma(\mathbf{D}).
$$ (3.8)
On account of conditions (3.3b, 3.3c), all elements of the matrix $\mathbf{D}$ are diagonal submatrices. The multidiagonal matrices $\frac{\partial F_{\mathrm{v}}}{\partial \mathbf{e}} - \frac{\partial \tilde{F}_{\mathrm{v}}}{\partial \mathbf{e}}$ and $\frac{\partial F_{\mathrm{v}}}{\partial \alpha}$ that are contained in the block $\vec{\mathbf{B}}$ are eliminated from the eigenvalue problem (3.8) due to the Schur determinant theorem (3.7) and the vector of zero blocks $\vec{\mathbf{0}}$.

<!-- Because of the diagonality of the blocks constituting $\mathbf{D}$, the remaining eigenvalue problem could be solved locally and in parallel for each grid point $i$:
$$ %\label{eq:jacspectrum2}
	\jmath(\mathbf{J}_\mathrm{ex,i}) = \{ 0 \} \cup \jmath(\mathbf{D}) \\
	= \{ 0 \} \cup \{ \jmath\left(\begin{bmatrix}
		  \frac{\partial F_{\mathrm{e},i}}{\partial \mathbf{e}_i}
		& \frac{\partial F_{\mathrm{e},i}}{\partial \alpha_i}
		& \frac{\partial F_{\mathrm{e},i}}{\partial \beta_i} \\[.7em]
		  \frac{\partial F_{\mathrm{a},i}}{\partial \mathbf{e}_i}
		& \frac{\partial F_{\mathrm{a},i}}{\partial \alpha_i} - \frac{\partial \tilde{F}_{\mathrm{a},i}}{\partial \alpha_i}
		& 0 \\[.7em]
		  \frac{\partial F_{\mathrm{b},i}}{\partial \mathbf{e}_i}
		& \frac{\partial F_{\mathrm{b},i}}{\partial \alpha_i}
		& \frac{\partial F_{\mathrm{b},i}}{\partial \beta_i} - \frac{\partial \tilde{F}_{\mathrm{b},i}}{\partial \beta_i}
	\end{bmatrix}\right), i = 1, \ldots, n \}.
$$
We use a local Rayleigh quotient iteration using the assembled adjugate of $\mathbf{J}_\mathrm{ex,i}$ to determine the spectral radius of the complete explicit Jacobian $\mathbf{J}_\mathrm{ex}$.
The diagonality properties listed above also mean that the spectrum of the explicit Jacobian $\mathbf{J}_\mathrm{ex}$ is independent of the distance associated with the spatial discretization, which fortunately means that the number of explicit stages is  to first order independent of grid refinement (but see [Section Y] for a more thorough analysis). -->

**3.2: Runge-Kutta extensions of the TR-BDF2 scheme**

We base our time integration scheme on the established TR-BDF2 scheme [references] with its explicit extension [Giraldo et al., 2013]. The TR-BDF2 scheme consists of a fractional trapezoidal (TR) step as a first stage, which is then completed with a second-order Backward Difference Formula (BDF2) stage. For an ODE $\partial y/\partial t = f(t, y)$ the scheme can be expressed as
$$
	y_{n+\gamma} - \left(\frac{\gamma}{2}\right) h_t f(t_n + \gamma h_t, y_{n+\gamma}) = y_n + \left(\frac{\gamma}{2}\right) h_t f(t_n, y_n),
$$ (3.9a)
$$
	y_{n+1} - \left(\frac{1-\gamma}{2-\gamma}\right) h_t f(t_n + h_t, y_{n+1}) = \left(\frac{1}{\gamma(2-\gamma)}\right) y_{n+\gamma} + \left(1 - \frac{1}{\gamma(2-\gamma)}\right) y_n.
$$ (3.9b)
Here, $h_t$ is the time step size, and $\gamma \in (0,1]$ is the fraction of the time step over which the first trapezoidal stage is computed, with a choice of $\gamma = 1$ representing a purely trapezoidal end member.
<!-- For a test equation $\partial y/\partial t = \lambda y$, $\lambda \in \mathbb{C}$, the TR-BDF2 scheme has the direct update form
\begin{subequations} \label{eq:trbdf2poly0} \begin{align}
	P_\mathrm{im}(h_t \lambda) y_{n+1} &= P_\mathrm{im,2}(h_t \lambda) P_\mathrm{im,1}(h_t \lambda) y_{n+1} = P_\mathrm{ex}(h_t \lambda) y_n, \\
	% \left(1 - \left(\frac{\gamma}{2}\right) h_t \lambda \right) \left(1 - \left(\frac{1-\gamma}{2-\gamma}\right) h_t \lambda \right) y_{n+1} &= \left(\frac{1}{\gamma(2-\gamma)}\right) \left(1 + \left(\frac{\gamma}{2}\right) h_t \lambda \right) y_n + \left(1 - \frac{1}{\gamma(2-\gamma)} \right) \left(1 - \left(\frac{\gamma}{2}\right) h_t \lambda \right) y_n.
	P_\mathrm{im,1}(z) &= \left(1 - \left(\frac{\gamma}{2}\right) z \right) \\
	P_\mathrm{im,2}(z) &= \left(1 - \left(\frac{1-\gamma}{2-\gamma}\right) z \right) \\
	P_\mathrm{im}(z) = P_\mathrm{im,2}(z) P_\mathrm{im,1}(z) &= 1 + \left( \frac{\gamma^2-2}{4-2\gamma} \right) z + \left( \frac{\gamma-\gamma^2}{4-2\gamma} \right) z^2 \\
	P_\mathrm{ex}(z) &= 1 + \left( \frac{1}{2-\gamma}-\frac{\gamma}{2}\right) z.
\end{align} \end{subequations}
We can see that the scheme is L-stable as long as $\gamma \neq 1$, since then the polynomial order of $P_\mathrm{im} > P_\mathrm{ex}$ and thus $\lim\limits_{h_t \lambda \to\pm\infty} P_\mathrm{ex}(h_t \lambda)/P_\mathrm{im}(h_t \lambda) = 0$ (i.e. for negative $\lambda$, the solution $y$ has a steady state on which it converges as the step size increases).

Plugging in the analytical solution $y_n = 1$, $y_{n+1} = \exp{h_t \lambda}$ and substituting the latter with its Taylor series approximation around $h_t \lambda = 0$, we obtain
\begin{subequations} \label{eq:trbdf2poly1} \begin{align}
	P_\mathrm{im}(h_t \lambda) \left(1 + h_t \lambda + \frac{1}{2}(h_t \lambda)^2 + \mathcal{O}(h_t^3) \right) - P_\mathrm{ex}(h_t \lambda) &= \mathcal{O}(h_t^3),
\end{align} \end{subequations}
which proves the second-order consistency of the combined scheme.

The free parameter $\gamma$ is chosen such that for the system of ODEs $\partial y / \partial t = J y$, the operators $P_\mathrm{im,2}(h_t J)$ and $P_\mathrm{im,1}(h_t J)$ are the same and are assembled only once per time step\footnote{We implement matrix-free algorithms and so the assembly cost is no issue for us.}. Thus,
\begin{align}
	\frac{\gamma}{2} = \frac{1-\gamma}{2-\gamma} \implies \gamma = 2 \pm \sqrt{2}.
\end{align}
The solution with minus sign, $2 - \sqrt{2} \approx 0.586$ not only lies in the desired interval $(0,1]$, but also yields a much lower coefficient of the $\mathcal{O}(h_t^3)$ truncation error ($\sim 0.04$ vs. $\sim 1.4$) and is therefore selected.

We continue the analysis of the TR-BDF2 scheme in the framework of a three-stage Runge-Kutta method, which is written for the autonomous ODE $\partial y/\partial t = f(y)$ as
\begin{subequations} \label{eq:rk2} \begin{align}
	w_1 &= y_n \\
	w_2 &= y_n + h_t\left[ \left(1-\frac{1}{\sqrt{2}}\right) f(w_1) + \left(1-\frac{1}{\sqrt{2}}\right) f(w_2) \right], \\
	w_3 &= y_n + h_t \left[ \left(\frac{1}{2\sqrt{2}}\right) f(w_1) + \left(\frac{1}{2\sqrt{2}}\right) f(w_2) + \left(1-\frac{1}{\sqrt{2}}\right) f(w_3) \right], \\
	y_{n+1} &= w_3.
\end{align} \end{subequations}
The first and last assignments are redundant, and the corresponding fields $w_1$ and $w_3$ do not need to be stored. We can however make use of the structure of \cref{eq:rk2} to create the derived method
\begin{subequations} \label{eq:rk3} \begin{align}
	w_1 &= y_n \\
	w_2 &= y_n + h_t\left[ \left(1-\frac{1}{\sqrt{2}}\right) f(w_1) + \left(1-\frac{1}{\sqrt{2}}\right) f(w_2) \right], \\
	w_3 &= y_n + h_t \left[ \left(\frac{1}{2\sqrt{2}}\right) f(w_1) + \left(\frac{1}{2\sqrt{2}}\right) f(w_2) + \left(1-\frac{1}{\sqrt{2}}\right) f(w_3) \right], \\
	y^\ast_{n+1} &= y_n + h_t \left[ b_1 f(w_1) + b_2 f(w_2) + b_3 f(w_3) \right],
\end{align} \end{subequations}
with coefficients $b^\mathrm{T} = [b_1, b_2, b_3]$, $\lVert b \rVert_1 = b_1 + b_2 + b_3 = 1$ that are to be chosen to yield a method that eliminates the third-order truncation error $\mathcal{O}(h_t^3)$, which does however come at the cost of L and A stability. Inserting again the test equation $\partial y/\partial t = \lambda y$ and following the procedure that led to \cref{eq:trbdf2poly1}, we now obtain
\begin{subequations} \label{eq:trbdf2poly2} \begin{align}
	P^\ast_\mathrm{im}(h_t \lambda) \left(1 + h_t \lambda + \frac{1}{2}(h_t \lambda)^2 + \frac{1}{6}(h_t \lambda)^3 + \mathcal{O}(h_t^4) \right) - P^\ast_\mathrm{ex}(h_t \lambda) &= g_1(b) (h_t \lambda) + g_2(b) (h_t \lambda)^2 + g_3(b) (h_t \lambda)^3 + \mathcal{O}(h_t^4),
\end{align} \end{subequations}
with the appropriate polynomials $P^\ast_\mathrm{im}$ and $P^\ast_\mathrm{ex}$, and coefficients $g(b)^\mathrm{T} = [g_1(b), g_2(b), g_3(b)]$, which are equated to zero to yield:
\begin{subequations} \label{eq:trbdf2poly3} \begin{align}
	P^\ast_\mathrm{im}(h_t \lambda) &\left(1 + h_t \lambda + \frac{1}{2}(h_t \lambda)^2 + \frac{1}{6}(h_t \lambda)^3 + \mathcal{O}(h_t^4) \right) - P^\ast_\mathrm{ex}(h_t \lambda) =  \mathcal{O}(h_t^4), \\
	b_1 &= \tfrac{1}{3} - \tfrac{1}{12}\sqrt{2}, \\
	b_2 &= \tfrac{1}{3} + \tfrac{1}{4} \sqrt{2}, \\
	b_3 &= \tfrac{1}{3} - \tfrac{1}{6} \sqrt{2}.
\end{align} \end{subequations}

Having computed a stable second-order accurate approximation $y_{n+1}$ to the solution $y(t_n + h_t)$, and a third-order accurate approximation $y^\ast_{n+1}$ to the same, we may estimate the $\mathcal{O}(h^3)$ truncation error $\epsilon_{n+1}$ as
\begin{align} \label{eq:trbdf2error1}
	\epsilon_{n+1} - \left(1-\frac{1}{\sqrt{2}}\right) h_t f(\epsilon_{n+1}) = y^\ast_{n+1} - y_{n+1}.
\end{align}
The implicit filter is due to [Shampine 1984, Hosea \& Shampine, 1996; in Bonaventura 2018] and gives the error measure $\epsilon^{n+1}$ the correct limit behavior for large $h_t$.

Finally, we seek a three-stage explicit Runge-Kutta scheme that complements the TR-BDF2 scheme by utilizing the same fractional step $\gamma h_t = (2 - \sqrt{2}) h_t$ and the same final stage coefficients $b^\mathrm{T} = [1/(2\sqrt{2}), 1/(2\sqrt{2}), 1 - 1/\sqrt{2}]$ (compare \cref{eq:rk2,eq:rk3}; see [Giraldo, 2013]). We thus seek the coefficient $a_{32}$ in
\begin{subequations} \label{eq:rk4} \begin{align}
	w_1 &= y_n \\
	w_2 &= y_n + h_t \left(2-\sqrt{2}\right) f(w_1), \\
	w_3 &= y_n + h_t \left[ (1 - a_{32}) f(w_1) + a_{32} f(w_2) \right], \\
	y^\ast_{n+1} &= y_n + h_t \left[ \left(\frac{1}{2\sqrt{2}}\right) f(w_1) + \left(\frac{1}{2\sqrt{2}}\right) f(w_2) + \left(1 - \frac{1}{\sqrt{2}}\right) f(w_3) \right],
\end{align} \end{subequations}
that maximizes the accuracy of the method. We again construct the update for the test equation $\partial y/\partial t = \lambda y$, $\lambda \in \mathbb{C}$,
\begin{align} \label{eq:expoly}
	y_{n+1} &= P(h_t \lambda) y_n, \\
	P(z) &= 1 + z + \frac{1}{2}z^2 + (3 - 2\sqrt{2}) a_{32} z^3 \\
	P(h_t \lambda) - \exp(h_t \lambda) &= \left((3 - 2\sqrt{2}) a_{32} - \frac{1}{6} \right) (h_t \lambda)^3 + \mathcal{O}(h_t^4) \\
	&= \mathcal{O}(h_t^4) \iff a_{32} = \tfrac{1}{2} + \tfrac{1}{3}\sqrt{2}.
\end{align}

\begin{figure}[htbp]
  \centering
  \label{fig:stabilityregions}\includegraphics[width=\columnwidth]{figures/stabilityregions}
  \caption{Stability regions $\lvert P(h_t \lambda) \rvert \leq 1$ in the complex plane $h_t \lambda \in \mathbb{C}$ derived from the scalar test equation $\partial y / \partial t = \lambda y$, of (\emph{a}) the second-order L-stable TR-BDF2 scheme, (\emph{b}) its third-order diagonally implicit Runge-Kutta extension, and (\emph{c}) its associated third-order explicit Runge-Kutta counterpart. In (\emph{c}), a circular region $D$ is plotted with origin $0$ and radius $\sqrt{3}$, the intersection of which with the left half plane is a proper subset of the stability region $\lvert P(h_t \lambda) \rvert \leq 1$ of the explicit Runge-Kutta scheme.}
  \label{fig:testfig}
\end{figure}

The stability regions of the three schemes derived here (second-order L-stable TR-BDF2, its third-order diagonally implicit Runge-Kutta extension, and its associated third-order explicit Runge-Kutta counterpart) are plotted in \cref{fig:stabilityregions}. In panel \emph{c} of the same figure, the region of stability of the explicit RK scheme is approximated with a circular region with radius $\sqrt{3}$ centered around the origin. The intersection of this region with the left half complex plane is a proper subset of the complete stability region, and allows us to determine an appropriate time step using only the the spectral radius of the explicit Jacobian $\mathbf{J}_\mathrm{ex}$ (see \cref{sect:imexsplit}):
\begin{align} \label{eq:timestep1}
	h_t = \sqrt{3}/\rho(\mathbf{J}_\mathrm{ex}).
\end{align}
Since the circular region intersects the imaginary axis at critical stability and passes close by the three complex zeroes of the stability polynomial, this choice gives both excellent conservation properties of the purely hyperbolic components, and excellent damping properties of the stiff dissipative components. We will refine the issue of time step selection in \cref{sect:errorcontrol}.

\subsection{TR-BDF2 in an Additive Runge-Kutta IMEX Scheme}

An idea first developed in [Giraldo, 2013], the additive IMEX split system \cref{eq:imexsystem} can be solved by a second-order Additive Runge-Kutta (ARK) discretization [references -- Giraldo2013 and references 1, 21, 26 therein], in which the third-order explicit scheme \cref{eq:rk4} is used for the explicit terms, the second-order, L-stable implicit scheme \ref{eq:rk2} is used for the implicit terms, and the second-order truncation error of the implicitly solved terms is estimated using \cref{eq:rk3,eq:trbdf2poly3,eq:trbdf2error1}.

The resulting Additive Runge-Kutta method is written as
\begin{subequations} \label{eq:imexark} \begin{align}
	w_i = y^n &+ h_t \sum\limits_{j=1}^{i-1} a_{ij} \left[ F(w_j) - \tilde{F}(w_j) \right] \nonumber \\
	          &+ h_t \sum\limits_{j=1}^{i} \tilde{a}_{ij} \tilde{F}(w_j), \qquad i = 1, \ldots, s \\
	y^{n+1} = y^n &+ h_t \sum\limits_{j=1}^{s} b_{j} F(w_j) \\
	(\mathbf{I} - \tilde{b}_{s} h_t \mathbf{J}_\mathrm{im}) \epsilon^{n+1} &= h_t \sum\limits_{j=1}^{s} (b^\ast_{j} - b_{j}) \tilde{F}(w_j), \label{eq:imexarkc}
\end{align} \end{subequations}
with $y^n = [ \vec{v}(t_n), \mathbf{e}(t_n), \alpha(t_n), \beta(t_n) ]^\mathrm{T}$, $F = [ F_\mathrm{v}, F_\mathrm{e}, F_\mathrm{a}, F_\mathrm{b} ]^\mathrm{T}$, $\tilde{F} = [ \tilde{F}_\mathrm{v}, \tilde{F}_\mathrm{e}, \tilde{F}_\mathrm{a}, \tilde{F}_\mathrm{b} ]^\mathrm{T}$, and coefficients $a$, $\tilde{a}$, $b$, $b^\ast$ given, based on the values derived in \cref{sect:trbdf2coeff}, by
\begin{subequations} \label{eq:coeff} \begin{align}
	a &= \begin{bmatrix}
		0                   & 0                & 0              \\[.7em]
		2-\sqrt{2}          & 0                & 0              \\[.7em]
		\frac{1}{2} - \frac{1}{3}\sqrt{2}  & \tfrac{1}{2} + \frac{1}{3}\sqrt{2} & 0              
	\end{bmatrix} \\
	\tilde{a} &= \begin{bmatrix}
		0                   & 0                & 0              \\[.7em]
		1 - \frac{1}{2}\sqrt{2}      & 1 - \frac{1}{2}\sqrt{2}   & 0              \\[.7em]
		\frac{1}{4}\sqrt{2}          & \frac{1}{4}\sqrt{2}       & 1 - \frac{1}{2}\sqrt{2} 
	\end{bmatrix} \\
	b &= \left[
		\tfrac{1}{4}\sqrt{2},\;\; \tfrac{1}{4}\sqrt{2},\;\; 1 - \tfrac{1}{2}\sqrt{2}
	\right]^\mathrm{T} \\
	b^\ast &= \left[
		\tfrac{1}{3} - \tfrac{1}{12}\sqrt{2},\;\; \tfrac{1}{3} + \tfrac{1}{4} \sqrt{2},\;\; \tfrac{1}{3} - \tfrac{1}{6} \sqrt{2}
	\right]^\mathrm{T}
\end{align} \end{subequations}

The overall truncation error of this ARK scheme is $\mathcal{O}(h_t^3)$ (it is thus accurate/consistent to second order), but the error measure $\epsilon_{n+1}$ only estimates the truncation error associated with the implicit components, since the explicit scheme already has the maximal third-order accuracy that can be expected of an explicit three-stage RK method. We do not however see this as a major drawback since the explicitly solved subsystem can be expected to be highly stiff, with the time step chosen by \cref{eq:timestep1} sufficiently affected by the unphysical peripheral components of the spectrum that it causes the physical components to be exceedingly accurately solved. Since the implicitly solved PDE components are not subject to a stability barrier, the same argument does not extend to them and their error must be monitored, and the time step further reduced if necessary.

\subsection{Error Control}\label{sect:errorcontrol}

Inspired by the multi-rate extension of TR-BDF2 of [Bonaventura et al., 2018], we design the following much simplified mechanism for controlling the error on the implicitly solved components of \cref{eq:imexsystem}. We estimate the implicit error using \cref{eq:imexarkc}. Then, given the absolute and relative tolerances $\tau_\mathrm{a}$ and $\tau_\mathrm{r}$, we compute the dimensionless error measure
\begin{align}
	\eta_i = \frac{\epsilon_{n+1,i}}{\tau_\mathrm{r} \lvert y_{n+1,i} \rvert + \tau_\mathrm{a}}.
\end{align}
We then compute the $L_\infty$ norm of $\eta$ over the vector blocks $[\vec{v}_{n+1},\mathbf{e}_{n+1}]$ (since these two will be solved together in block reduced form), $[\alpha_{n+1}]$, and $[\beta_{n+1}]$. We denote these maximum dimensionless errors by $\eta_\mathrm{v}$, $\eta_\mathrm{a}$, and $\eta_\mathrm{b}$, respectively.

If any of these measures exceed a value of one, a new time step is computed by [Bonaventura 2018 and references [14], [19] therein]
\begin{align}
	h_t^\ast = \nu h_t \mathrm{max}(\eta_\mathrm{v}, \eta_\mathrm{a}, \eta_\mathrm{b})^{-1/3},
\end{align}
with $\nu \in (0, 1)$ a safety coefficient and the power $1/3$ specific to a method with second-order accuracy.

The time step is then redone with $h_t^\ast$ and this process is repeated until $\mathrm{max}(\eta_\mathrm{v}, \eta_\mathrm{a}, \eta_\mathrm{b}) \leq 1$. When this condition is satisfied, there ought to be some lingering memory of the failure of the explicit step, and the time step selection \cref{eq:timestep1} is best amended with the recursion
\begin{align} \label{eq:timestep2}
	h_t^{(n+1)} = \mathrm{min}( \sqrt{3}/\rho(\mathbf{J}_\mathrm{ex}), (1/\nu) h_t^{(n)} ),
\end{align}
with the same safety coefficient $\nu \in (0, 1)$.

At each stage of time step refinement the initial guess of the new solution may be significantly improved by exploiting the globally $C_1$-continuous cubic Hermite interpolation that is given in [Section 5 of Bonaventura et al., 2018], which requires only the known solution stages in its evaluation.

[It might actually be possible to save some refinements by solving a constraint optimization problem using the cubic Hermite interpolation of the solution over the time step, optimizing the non-dimensional error measure $\eta$ subject to the constraint $\eta \geq \nu$.] -->