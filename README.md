# Rayleigh Benard Convection

## Governing Equations

Incompressible Navier-Stokes equations with body force $- g\ \widehat{k}$ :

$$\mathbf{\nabla}.\mathbf{u} = 0$$ (1)

$$\rho\frac{\partial\mathbf{u}}{\partial t} + \rho\mathbf{{u}}.\mathbf{\nabla}\mathbf{u} = - \mathbf{\nabla}p - \rho g\widehat{k} + \mu\mathbf{\nabla}^{2}\mathbf{u}$$ (2)

where, $\mathbf{u} = (u,v,w)$ are the velocities in $x,\ y$ and $z$ directions respectively, $p$ is the pressure, $\rho$ is the density and $\mu$ is the dynamic viscosity.

## Boussinesq Approximation

This approximation is based on the fact that variation of density induced by the variation of temperature can be everywhere except in the buoyancy term.

Say, the density is expressed as,

$$\rho = \rho_{c}\left( 1 - \alpha\left( T - T_{c} \right) \right)$$ (3)

where, $\rho_{c}$ and $T_{c}$ are the density and the temperature of the fluid at the cold plate respectively.

Substituting $\rho$ in the momentum equation, we get,

$$\rho_{c}\left( 1 - \alpha\left( T - T_{c} \right) \right)\left( \frac{\partial\mathbf{u}}{\partial t} + \mathbf{u}.\mathbf{\nabla}\mathbf{u} \right) = - \mathbf{\nabla}p - \rho_{c}\left( 1 - \alpha\left( T - T_{c} \right) \right)g\widehat{k} + \mu\mathbf{\nabla}^{2}\mathbf{u}$$ (4)

Now, applying the Boussinesq approximation, we neglect the dependence of density in the inertia terms.

$$\rho_{c}\left( \frac{\partial\mathbf{u}}{\partial t} + \mathbf{u}.\mathbf{\nabla}\mathbf{u} \right) = - \mathbf{\nabla}p - \rho_{c}\left( 1 - \alpha\left( T - T_{c} \right) \right)g\widehat{k} + \mu\mathbf{\nabla}^{2}\mathbf{u}$$
$$\rho_{c}\left( \frac{\partial\mathbf{u}}{\partial t} + \mathbf{u}.\mathbf{\nabla}\mathbf{u} \right)\mathbf{= \ } - \left( \mathbf{\nabla}p + \rho_{c}g\widehat{k} \right) + \rho_{c}\alpha\left( T - T_{c} \right)g\widehat{k} + \mu\mathbf{\nabla}^{2}\mathbf{u}$$ (5)

Define:

$$p^{'} = p + \rho_{c}gz$$ (6)

Therefore,

$$\mathbf{\nabla}p^{'} = \ \mathbf{\nabla}p + \rho_{c}g\widehat{k}$$ (7)

Substituting (7) in (5), we get,

$$\rho_{c}\left( \frac{\partial\mathbf{u}}{\partial t} + \mathbf{u}.\mathbf{\nabla}\mathbf{u} \right)\mathbf{= \ } - \mathbf{\nabla}p^{'} + \rho_{c}\alpha\left( T - T_{c} \right)g\widehat{k} + \mu\mathbf{\nabla}^{2}\mathbf{u}$$ (8)

Since, the number of variables is more than the number of
equations, we need energy equation to find $T$.

$$\frac{\partial T}{\partial t} + \mathbf{u}.\mathbf{\nabla}T = \kappa\mathbf{\nabla}^{2}T$$ (9)

