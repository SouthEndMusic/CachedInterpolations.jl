# Inverting the integral

We are interested in inverting $U(t)$ as defined above. Note that $U(t)$ is invertible if $u(t)$ is positive for all $t > t_1$. If we define

```math
\begin{equation}
    U_J = \int_{t_1}^{\tilde{t}_J}u(\tau)\text{d}\tau = \sum_{j = 2}^{J} \int_{\tilde{t}_{j-1}}^{\tilde{t}_j} u(\tau)\text{d}\tau,
\end{equation}
```

then solving $U(t) = V$ for $t$ where $V \in [U_{J-1}, U_J]$ amounts to solving

```math
\begin{equation}
    \int_{\tilde{t}_{J-1}}^t u(\tau)\text{d}\tau = V - U_{J-1}.
\end{equation}
```

For linear sections this yields a quadratic equation in $t$ with solution

```math
\begin{equation}
    t = t_I + \left[-\frac{u_I}{\Delta u_{I+1}} + \text{sign}\left(\frac{u_{I+1}}{u_{I+1}}\right)\sqrt{\left(\frac{u_I}{\Delta u_{I+1}}\right)^2 +\lambda\left(\frac{u_I}{\Delta u_{I+1}} + \frac{\lambda}{4}\right) +2\frac{V - U_{J-1}}{\Delta t_{I+1}\Delta u_{I+1}}}\right]\Delta t_{I+1}.
\end{equation}
```

 For spline sections this leads to a quartic equation in $s$:

```math
\begin{equation}
    \begin{aligned}
        3(\Delta t_{I+1} - \Delta t_I)(\Delta u_{I+1} - \Delta u_I)s^4 + \\
        4\Delta t_I (\Delta u_{I+1} - \Delta u_I) s^3 + \\
        12(\Delta t_{I+1} - \Delta t_I)\left(\Delta u_I + \frac{u_{I - \frac{\lambda}{2}}}{\lambda}\right)s^2 + \\
        24 \Delta u_I \left(\Delta u_I + \frac{u_{I - \frac{\lambda}{2}}}{\lambda}\right) s + \\
        - \frac{24}{\lambda^2}(V - U_{J-1}) = 0
    \end{aligned}
\end{equation}
```
 
 
This quartic equation can be solved with the [quartic formula](https://en.wikipedia.org/wiki/Quartic_function#General_formula_for_roots).

[Quartic equation at Wolfram Alpha](https://www.wolframalpha.com/input?i=integrate+%28lambda*%28t_2+-+t_1%29*s+%2B+lambda*t_1%29*%28lambda%2F2+*+%28u_2+-+u_1%29*s%5E2+%2B+lambda*u_1+%2B+u_3%29+ds+from+0+to+S)
