This simple program is aimed to apply the finite distance approximation to the solution of the diffusion equation, useful for the solution of options' value.
## Finite distance approxiamtion
Suppose that we want to solve the diffusion equation ($u$ here representing for example the temperature in a long, thin and uniform bar of metal as function of space and time)

$$\frac{∂ u}{∂\tau}=\frac{∂^2 u}{∂ x^2},$$

with initial conditions 

$$u(x,0)=u_0(x)\hspace{0.5cm}u(x,\tau)\sim u_{-\infty}(x,\tau)\hspace{0.5cm}u(x,\tau)\sim u_{\infty}(x,\tau)\hspace{0.5cm}as\hspace{0.3cm}x\rightarrow\pm\infty.$$

Using the finite distance-approximation, one can show that the first equation becomes

$$u_{n}^{m+1}=\alpha u_{n+1}^{m}+(1-2\alpha)u_{n}^{m}+\alpha u^{m}_{n-1}\hspace{0.5cm}where\hspace{0.5cm}\alpha\equiv\frac{δ \tau}{(δ x)^2}\hspace{0.5cm}and\hspace{0.5cm}u_{n}^{m}(x,\tau)\equiv u(nx,m\tau).$$

In this file there is a program which iteratively solve this equation.
