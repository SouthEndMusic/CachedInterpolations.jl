# Integrating

## Complete intervals

We are interested in integrating the smoothed interpolation from the start to some $t > t_1$. To compute this efficiently we need to know the integral of the interpolation over the various intervals. More precisely:

- For the linear sections we obtain

```math
\begin{equation}
    \int_{t_{i-1+\frac{\lambda}{2}}}^{t_{i-\frac{\lambda}{2}}} 
         u_{i-1} + \frac{\Delta u_i}{\Delta t_i}(\tau - t_{i-1})
    \text{d}\tau
    =
    (1 - \lambda)\Delta t_i 
    \left[
        u_{i-1}
        +
        \frac{1}{2}(1 - \lambda) \Delta u_i
    \right].
\end{equation}
```

- For the spline sections we obtain

```math
\begin{equation}
    \begin{aligned}
        \int_{t_{i-\frac{\lambda}{2}}}^{t_{i + \frac{\lambda}{2}}} u_i(\tau)\text{d}\tau
        &=&
        \int_0^1 T'_i(s)u_i(T_i(s))\text{d}s \\ 
        &=&
        \int_0^1
            \left[\lambda(\Delta t_{i+1} - \Delta t_i)s + \lambda\Delta t_i\right]
            \left[\frac{\lambda}{2}(\Delta u_{i+1} - \Delta u_i)s^2 + \lambda \Delta u_i s + u_{i - \frac{\lambda}{2}}\right]
        \text{d}s \\
        &=& \frac{\lambda^2}{24}
        \left[
            \Delta t_i \left(-3\Delta u_i + \Delta u_{i+1}\right) +
            \Delta t_{i+1} \left(-\Delta u_i + 3 \Delta u_{i+1}\right)
        \right]
        +
        \frac{\lambda}{2}(\Delta t_i + \Delta t_{i+1})u_i.
    \end{aligned}
\end{equation}
```

## Incomplete intervals

We now define the new set of points $(\tilde{\mathbf{p}}_j)_{j=1}^{2n}$ given by all the $\mathbf{p}_{i - \frac{\lambda}{2}}, \mathbf{p}_{i+ \frac{\lambda}{2}}$ and the original boundary points, sorted by $t$. Then for $t \in \left[\tilde{t}_{J-1}, \tilde{t}_J\right]$ we have

```math
\begin{equation}
    \begin{aligned}
        U(t) = \int_{t_1}^t u(\tau)\text{d}\tau = \sum_{j = 2}^{J-1} \int_{\tilde{t}_{j-1}}^{\tilde{t}_j} u(\tau)\text{d}\tau + \int_{\tilde{t}_{J-1}}^t u(\tau)\text{d}\tau,
    \end{aligned}
\end{equation}
```

Where the summed integrals are given by the values above. For the last integral:

- If $J$ is odd then the last integral is of a linear section:

```math
\begin{equation}
    \int_{\tilde{t}_{J-1}}^t u(\tau)\text{d}\tau = \left((t-t_I) - \frac{\lambda}{2}\Delta t_{I+1}\right)u_I
    + 
    \frac{1}{2}\frac{\Delta u_{I+1}}{\Delta t_{I+1}}\left[(t-t_I)^2 - \frac{\lambda^2}{4}\Delta t_{I+1}^2\right]
\end{equation}
```

where $I = \frac{J-1}{2}$.

- If $J$ is even the last integral is of a spline section:

```math
\begin{equation}
    \begin{aligned}
        \int_{\tilde{t}_{J-1}}^t u_I(\tau)\text{d}\tau &=& 
        \int_0^{S_I(t)} T'_I(s)u_I\left(T_I(s)\right)\text{d}s \\
        &=&
        \int_0^{S_I(t)}
            \left[\lambda(\Delta t_{I+1} - \Delta t_I)s + \lambda\Delta t_I\right]
            \left[\frac{\lambda}{2}(\Delta u_{I+1} - \Delta u_I)s^2 + \lambda \Delta u_I s + u_{I - \frac{\lambda}{2}}\right]
        \text{d}s
    \end{aligned}
\end{equation}
```

where $I = \frac{J}{2}$.
