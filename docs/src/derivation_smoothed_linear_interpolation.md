# Derivation of smoothed linear interpolation

## Linear interpolation

Given is the set of points $(\mathbf{p}_i)_{i=1}^n$ in $\mathbb{R}^2$, where we write $t$ and $u$ for the respective coordinates. We assume that the $t_i$ are strictly increasing. Linear interpolation of these points is then simply given by

```math
\begin{equation}
    u|_{[t_{i-1}, t_i]}(t) = \frac{u_i-u_{i-1}}{t_i-t_{i-1}}(t-t_{i-1}).
\end{equation}
```

Note that there is a discontinuity in derivative of the function $u$ at each $t_i$ for $i = 2, \ldots, n-1$.

## Smoothed spline corners

To get rid of the discontinuities mentioned in the previous section, we take out a section of the interpolation around each discontinuity and replace it with a spline curve.

### New points

For the construction of the smoothing we consider the consecutive points

```math
    \mathbf{p}_{i-1}, \mathbf{p}_{i}, \mathbf{p}_{i+1}.
```

Now we disregard $\mathbf{p}_i$, and introduce 2 new points

```math
\begin{equation}
    \begin{aligned}
        \mathbf{p}_{i-\frac{\lambda}{2}} =&\; \mathbf{p}_i - \frac{\lambda}{2}\Delta\mathbf{p}_i \\
        \mathbf{p}_{i+\frac{\lambda}{2}} =&\; \mathbf{p}_i + \frac{\lambda}{2}\Delta\mathbf{p}_{i+1}
    \end{aligned}
\end{equation}
```

where $\Delta\mathbf{p}_i = \mathbf{p}_i - \mathbf{p}_{i-1}$ and $\lambda \in [0,1]$. We will connect these points with a spline curve, and so $\lambda$ determines the size of the interval around $\mathbf{p}_i$ that is replaced by the spline curve.

### Deriving the spline curve

We want to connect $\mathbf{p}_{i-\frac{\lambda}{2}}$ and $\mathbf{p}_{i+\frac{\lambda}{2}}$ with a smooth parametric curve 
```math
\begin{equation}
    \mathbf{C}_i : [0,1] \rightarrow \mathbb{R}^2
\end{equation}
```

such that:

- The connection can be expressed as 

```math
\begin{equation}
    u_i : \left[t_{i - \frac{\lambda}{2}}, t_{i + \frac{\lambda}{2}}\right] \rightarrow \mathbb{R}, 
\end{equation}
```

i.e. the $t$ component of the curve must be invertible.

- The connection is continuous, i.e.

```math
\begin{equation}
    \mathbf{C}_i(0) = \mathbf{p}_{i-\frac{\lambda}{2}}, \quad \mathbf{C}_i(1) = \mathbf{p}_{i+\frac{\lambda}{2}}.
\end{equation}
```

- The derivative of the connection is continuous, i.e.
```math
\begin{equation}
    \mathbf{C}'_i(0) \propto \Delta\mathbf{p}_i, \quad \mathbf{C}'_i(1) \propto \Delta\mathbf{p}_{i+1}.
\end{equation}
```

We can achieve this by repeated interpolation. The first interpolations are 

```math
\begin{equation}
    \begin{aligned}
        \mathbf{C}_{i-\frac{\lambda}{2}}(s) = (1-s)\mathbf{p}_{i-\frac{\lambda}{2}} + s\mathbf{p}_i \\
        \mathbf{C}_{i+\frac{\lambda}{2}}(s) = (1-s)\mathbf{p}_i + s\mathbf{p}_{i+\frac{\lambda}{2}}
    \end{aligned}
\end{equation}
```

and combining these yields

```math
\begin{equation}
    \begin{aligned}
        \mathbf{C}_i(s) =&\; (1-s)\mathbf{C}_{i-\frac{\lambda}{2}}(s) + s\mathbf{C}_{i+\frac{\lambda}{2}}(s) \\
        =&\; (1-s)^2\mathbf{p}_{i-\frac{\lambda}{2}} + 2s(1-s)\mathbf{p}_i + s^2\mathbf{p}_{i+\frac{\lambda}{2}} \\
        =&\; \frac{\lambda}{2}(\Delta \mathbf{p}_{i+1} - \Delta \mathbf{p}_i)s^2 + \lambda \Delta \mathbf{p}_i s + \mathbf{p}_{i-\frac{\lambda}{2}}
    \end{aligned}
\end{equation}
```

Note that the second formulation tells us that $C_i$ is a convex combination of $\mathbf{p}_{i-\frac{\lambda}{1}}, \mathbf{p}_i, \mathbf{p}_{i + \frac{\lambda}{2}}$ for all $s \in [0,1]$ and thus always is in the convex hull of these points.

### Writing spline curve as a function $u(t)$

To write the spline curve as a function $u(t)$, we first need to obtain $s$ from $t$:

```math
\begin{equation}
    T_i(s) = \frac{1}{2}\lambda(\Delta t_{i+1} - \Delta t_i)s^2 + \lambda\Delta t_i s + t_{i-\frac{\lambda}{2}} = t.
\end{equation}
```

This yields

```math
\begin{equation}
    S_i(t) = \frac{
            -\lambda \Delta t_i + \sqrt{\lambda^2\Delta t_i^2 + 2\lambda (\Delta t_{i+1} - \Delta t_i)\left(t - t_{i-\frac{\lambda}{2}}\right)}
        }{
            \lambda (\Delta t_{i+1} - \Delta t_i)
        },
\end{equation}
```

or, in the degenerate case that $\Delta t_{i+1} - \Delta t_i = 0$ (i.e. the 3 points are equally spaced),

```math
\begin{equation}
    S_i(t) = \frac{1}{\lambda}\frac{t - t_{i - \frac{\lambda}{2}}}{\Delta t_i}.
\end{equation}
```

Note that $\Delta t_i \ne 0$ by the assumption that the $t_n$ are strictly increasing.

We conclude:

```math
\begin{equation}
    u_i(t) = \frac{\lambda}{2}(\Delta u_{i+1} - \Delta u_i)S_i(t)^2 + \lambda \Delta u_i S_i(t) + u_{i - \frac{\lambda}{2}}.
 \end{equation}
```

## Extrapolation

## Evaluation