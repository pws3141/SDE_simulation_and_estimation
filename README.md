# TODO: README is old, some information is incorrect


# Ornstein-Uhlenbeck Simulation and Parameter Estimation
This repo contains R functions to do the following stochastic differential equation
simulation and parameter estimation for either $X_t$ or $Y_t$ , $t \in [0, T]$, 
given by the equations
> dX_t = - \gamma X_t dt + \sigma dB_t
> Y_t = X_t + m(t)
where $m(t)$ is a continuous function given by
> m(t) = \sum_{i = 1}^I m_i(t) \mathbb{I}_{t \in B_i}
where $m_i(t), i = 1, 2, \ldots, I$ are quadratic functions and 
$\{B_i, i = 1, 2, \ldots, I\}$ is a partition of the interval [0, T), closed from
the left.

## Part 1: Splines
The splines $m_t(t), i = 1, 2, \ldots, I$ are quadratic functions found using the
following three constraints:
        1) $\int_{B_j} m_j(t) dt = \mu_j$ for some specified $\mu_j$;
        2) $m_{j-1}(t) = m_j(t)$ for $t = \min\{t: t \in B_j\}$;
        3) $\frac{d}{dt} m_{j-1}(t) = \frac{d}{dt} m_j(t)$ for $t = \min\{t: t \in B_j\}$.

Note here that constraint (1) uses the mean over each partition $B_i$ and not just a
singular point to caluculate the splines.

The splines can be estimated from data by taking the arithmetic mean of all points that belong
to each interval $B_i, i = 1, 2, \ldots, I$, which can then be used in constraint
(1).

## Part 2: Simulation
Simulation is done using Euler-Maruyama method. Note that if we want to simulate
$Y_t$ we have to use both $m(t)$ and $\frac{d}{dt} m(t)$ as
> dY_t = (-\gamma (X_t - m(t)) + \frac{d}{dt} m(t)) dt + \sigma dB_t.

## Part 3: Parameter Estimation using Maximum Likelihood
Maximum likelihood is performed on a sparsely observed chain
$\{X_i, i = 1, 2, \ldots, n}$, where the obsrervations are at regularly
spaced time intervals $\delta t$. If we observe $\{Y_t, t = 0, \delta t, \ldots, T\}$ 
then we can remove the splines $m(t)$ before estimating the parameters.

As $X_t$ is assumed to  be an Ornstein-Uhlenbeck process we therefore know the
form of the transition density $p(t; x, y)$, and the log-likelihood is given by
> l(X; gamma, sigma^2) = \sum_{i = 1}^{n - 1} p(\delta t; X_i, X_{i + 1})

The maximum likelihood estimators $\hat \gamma$ and $\hat \sigma^2$ are found by
maximising the log-likelihood using the 'optim' function in R. As 
$\hat \sigma^2 = g(\hat \gamma)$ we only need to optimise over possible $\gamma$
values.

### Remark
Note that the accuracy of the parameter estimation method increases in the limit as 
$n \rightarrow \infty$ with $\delta t$ constant. If we simulate a chain using Step
(2) and then want to perform MLE as described above, the chain should first be
'thinned' by taking every $N$ observations such that the new 'sparse' chain is given
by $X_0, X_{N \delta t}, X_{2N \delta t}, \ldots, X_T$.
