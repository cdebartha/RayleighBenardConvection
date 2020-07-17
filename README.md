# Rayleigh Benard Convection

## Governing Equations

Incompressible Navier-Stokes equations with body force <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/930600685a1f78702b6e7cc965660e2d.svg?invert_in_darkmode" align=middle width=35.799139199999985pt height=32.42016360000002pt/> :

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/719f23136027ff8435fc8451db709992.svg?invert_in_darkmode" align=middle width=58.90396379999999pt height=11.232861749999998pt/></p>

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/e2524fbcd54b0d0bc01df9da60a4d79a.svg?invert_in_darkmode" align=middle width=264.78570525pt height=33.81208709999999pt/></p>

where, <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/3d96202da71b051caff14815ac3573bd.svg?invert_in_darkmode" align=middle width=89.99601764999998pt height=24.65753399999998pt/> are the velocities in <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/86ad5e70a1497f100df341d42edabbdd.svg?invert_in_darkmode" align=middle width=30.82954049999999pt height=14.15524440000002pt/> and <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/f93ce33e511096ed626b4719d50f17d2.svg?invert_in_darkmode" align=middle width=8.367621899999993pt height=14.15524440000002pt/> directions respectively, <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/2ec6e630f199f589a2402fdf3e0289d5.svg?invert_in_darkmode" align=middle width=8.270567249999992pt height=14.15524440000002pt/> is the pressure, <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?invert_in_darkmode" align=middle width=8.49888434999999pt height=14.15524440000002pt/> is the density and <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/07617f9d8fe48b4a7b3f523d6730eef0.svg?invert_in_darkmode" align=middle width=9.90492359999999pt height=14.15524440000002pt/> is the dynamic viscosity.

## Boussinesq Approximation

This approximation is based on the fact that variation of density induced by the variation of temperature can be everywhere except in the buoyancy term.

Say, the density is expressed as,

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/dcd7ec0bd28b17c9fa06026d90e628b6.svg?invert_in_darkmode" align=middle width=163.83224879999997pt height=16.438356pt/></p>

where, <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/4f649235e80d795f0a153fc88d1e7d44.svg?invert_in_darkmode" align=middle width=14.373536099999988pt height=14.15524440000002pt/> and <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/f3087bed09ac6d957af8ecfe91212479.svg?invert_in_darkmode" align=middle width=15.480837899999988pt height=22.465723500000017pt/> are the density and the temperature of the fluid at the cold plate respectively.

Substituting <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?invert_in_darkmode" align=middle width=8.49888434999999pt height=14.15524440000002pt/> in the momentum equation, we get,

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/de1fee5091fd8a35af2e25cf18e6618b.svg?invert_in_darkmode" align=middle width=535.8008238pt height=39.452455349999994pt/></p>

Now, applying the Boussinesq approximation, we neglect the dependence of density in the inertia terms.

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/2fc49f8455f21ed28b23b830766e1fcb.svg?invert_in_darkmode" align=middle width=417.58055955pt height=39.452455349999994pt/></p>
<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/62a04feb41e0b6fc9b7dca0ac4626f97.svg?invert_in_darkmode" align=middle width=452.56479015pt height=39.452455349999994pt/></p>

Define:

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/5f4c1b3e7092df86604873fe590e1e43.svg?invert_in_darkmode" align=middle width=95.81168355pt height=19.5105966pt/></p>

Therefore,

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/8d922d8b36f7361fa3e46966a3b41efb.svg?invert_in_darkmode" align=middle width=129.4247526pt height=19.5105966pt/></p>

Substituting, we get,

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/16dad1c6b6914e6caa2bc032d29281ec.svg?invert_in_darkmode" align=middle width=385.40603805pt height=39.452455349999994pt/></p>

Since, the number of variables is more than the number of equations, we need energy equation to find <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/2f118ee06d05f3c2d98361d9c30e38ce.svg?invert_in_darkmode" align=middle width=11.889314249999991pt height=22.465723500000017pt/>.

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/3a816856923535f3eae73abaeda2e21e.svg?invert_in_darkmode" align=middle width=148.6011054pt height=33.81208709999999pt/></p>

## Non-dimensionalization

Define:

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/97047b975435d9e8c16095306833aec0.svg?invert_in_darkmode" align=middle width=203.4271569pt height=29.47417935pt/></p>

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/9dd9b22ad921795e6f392c15803d61c0.svg?invert_in_darkmode" align=middle width=259.35192524999997pt height=33.58376834999999pt/></p>

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/7a4b58162ee200c1816ae78bb499b168.svg?invert_in_darkmode" align=middle width=323.80429125pt height=44.59463415pt/></p>

where, <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/78ec2b7008296ce0561cf83393cb746d.svg?invert_in_darkmode" align=middle width=14.06623184999999pt height=22.465723500000017pt/> is the distance between the cold plate and the hot plate, <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/33c1b636b1f08477704e3ff994d985b1.svg?invert_in_darkmode" align=middle width=17.30221349999999pt height=22.465723500000017pt/> is the temperature of the fluid at the hot plate and <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/b49211c7e49541e500c32b4d56d354dc.svg?invert_in_darkmode" align=middle width=9.16670204999999pt height=14.15524440000002pt/> is the kinematic viscosity of the fluid.

## Non-dimensionalized equations

Ignoring <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/239f0fec4e07802db7bd8c11974f65e8.svg?invert_in_darkmode" align=middle width=8.219209349999991pt height=15.296829900000011pt/> from the non-dimensionalized variables, we get,

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/3257fa6f0842cb89d9578499c10d891f.svg?invert_in_darkmode" align=middle width=150.33260055pt height=37.0084374pt/></p>

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/8f154fd495df3c7164c8ba97cbcff187.svg?invert_in_darkmode" align=middle width=414.9561999pt height=38.973783749999996pt/></p>

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/5885ee08be4141976074a59260e6d2b8.svg?invert_in_darkmode" align=middle width=410.10040004999996pt height=38.973783749999996pt/></p>

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/1ccb50ff2d2434dc07cf8a33b8d3409b.svg?invert_in_darkmode" align=middle width=523.17621345pt height=40.11819404999999pt/></p>

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/8cc8cae8fefeaf5fbf35fd9dc4269319.svg?invert_in_darkmode" align=middle width=531.2389400999999pt height=40.11819404999999pt/></p>

where,

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/c96b750a458ab4a9c7bac2aa76ea6f0b.svg?invert_in_darkmode" align=middle width=158.41814175pt height=35.77743345pt/></p>
<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/33a17206b72a0e2f8d1038125817762f.svg?invert_in_darkmode" align=middle width=54.071069249999994pt height=29.47417935pt/></p>

## Computational domain

![Computational domain](images/domain.png)

## Computational grid

For testing the code, an initial grid of 64 × 64 × 32 is used before implementing the code on the main grid 128 × 128 × 64.

![grid64](images/grid64.png)
![grid128](images/grid128.png)

For discretizing the governing equations, staggered grid system is used where pressure and temperature are stored in the cell center and velocities are stored in cell faces.

    Nx = 64
    Ny = 64
    Nz = 32
    Lx = 8.0
    Ly = 8.0
    Lz = 1.0
    hx = Lx/(Nx-1)
    hy = Ly/(Ny-1)
    hz = Lz/(Nz-1)
    #Generate grid for cell edges
    xe = np.linspace(0,Lx,Nx)
    ye = np.linspace(0,Ly,Ny)
    ze = np.linspace(0,Lz,Nz)
    x, y, z = np.meshgrid(xe, ye, ze, indexing='ij')
    #Generate grid for cell centers
    xc = np.linspace(-hx/2,Lx+hx/2,Nx+1)
    yc = np.linspace(-hy/2,Ly+hy/2,Ny+1)
    zc = np.linspace(-hz/2,Lz+hz/2,Nz+1)

## Initial Conditions

All the variables <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/>, <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/>, <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/31fae8b8b78ebe01cbfbe2fe53832624.svg?invert_in_darkmode" align=middle width=12.210846449999991pt height=14.15524440000002pt/>, <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/2ec6e630f199f589a2402fdf3e0289d5.svg?invert_in_darkmode" align=middle width=8.270567249999992pt height=14.15524440000002pt/> and <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/2f118ee06d05f3c2d98361d9c30e38ce.svg?invert_in_darkmode" align=middle width=11.889314249999991pt height=22.465723500000017pt/> are given an initial value of 0 at all the grid points.

    u = np.zeros((Nx,Ny+1,Nz+1))
    v =  np.zeros((Nx+1,Ny,Nz+1))
    w = np.zeros((Nx+1,Ny+1,Nz))
    p = np.zeros((Nx+1,Ny+1,Nz+1))
    T = np.zeros((Nx+1,Ny+1,Nz+1))

where, `Nx`, `Ny`, `Nz` are the number of grid points in x, y and z directions respectively.

## Boundary Conditions

For x and y directions, periodic boundary conditions are applied.

    # Periodic in x
    for j in range(Ny+1):
        for k in range(Nz+1):
            u[0,j,k] = u[Nx-2,j,k]
            u[Nx-1,j,k] = u[1,j,k]
            p[0,j,k] = p[Nx-1,j,k]
            p[Nx,j,k] = p[1,j,k]
            T[0,j,k] = T[Nx-1,j,k]
            T[Nx,j,k] = T[1,j,k]
            if j != Ny:
                v[0,j,k] = v[Nx-1,j,k]
                v[Nx,j,k] = v[1,j,k]
            if k != Nz:
                w[0,j,k] = w[Nx-1,j,k]
                w[Nx,j,k] = w[1,j,k]
    # Periodic in y
    for i in range(Nx+1):
        for k in range(Nz+1):
            v[i,0,k] = v[i,Ny-2,k]
            v[i,Ny-1,k] = v[i,1,k]
            p[i,0,k] = p[i,Ny-1,k]
            p[i,Ny,k] = p[i,1,k]
            T[i,0,k] = T[i,Ny-1,k]
            T[i,Ny,k] = T[i,1,k]
            if i != Nx:
                u[i,0,k] = u[i,Ny-1,k]
                u[i,Ny,k] = u[i,1,k]
            if k != Nz:
                w[i,0,k] = w[i,Ny-1,k]
                w[i,Ny,k] = w[i,1,k]

For z-direction, no-slip boundary condition is applied for velocities and pressure and fixed boundary condition is applied for temperature.

    # For velocities, no-slip at z-min and z-max
    # For temperature, all-fixed at z-min and zmax
    for i in range(Nx+1):
        for j in range(Ny+1):
            w[i,j,0] = 0
            w[i,j,Nz-1] = 0
            p[i,j,0] = p[i,j,1]
            p[i,j,Nz] = p[i,j,Nz-1]
            T[i,j,0] = 2 - T[i,j,1] # T(zmin) = 1
            T[i,j,Nz] = - T[i,j,Nz-1] # T(zmax) = 0
            if i != Nx:
                u[i,j,0] = - u[i,j,1]
                u[i,j,Nz] = - u[i,j,Nz-1]
            if j != Ny:
                v[i,j,0] = - v[i,j,1]
                v[i,j,Nz] = - v[i,j,Nz-1]

## Formulation of the discretized equations

Mixed RK3-CN time-stepping technique is used. Explicit terms are treated using RK3 Williamson scheme and the implicit terms are computed using Crank-Nicholson scheme. Since the grid size is significantly less in the z-direction, the double derivative terms in z-direction are treated as implicit terms.

## Calculation of the explicit terms

For u-momentum and v-momentum equations, explicit terms consist of the advection terms and the diffusion terms.

    Eu = np.zeros((Nx,Ny+1,Nz+1))
    Ev = np.zeros((Nx+1,Ny,Nz+1))
    Ew = np.zeros((Nx+1,Ny+1,Nz))
    ET = np.zeros((Nx+1,Ny+1,Nz+1))
    
    for k in range(1,Nz):
        for j in range(1,Ny):
            for i in range(1,Nx-1):
                # U2x
                uavg2 = 0.5*(u[i+1,j,k]+u[i,j,k])
                uavg1 = 0.5*(u[i,j,k]+u[i-1,j,k])
                
                U2x = ((uavg2)**2-(uavg1)**2)/(xc[i+1]-xc[i])
                
                #UVy
                vavg2 = 0.5*(v[i,j,k]+v[i+1,j,k])
                vavg1 = 0.5*(v[i,j-1,k]+v[i+1,j-1,k])
                uavg2 = 0.5*(u[i,j,k]+u[i,j+1,k])
                uavg1 = 0.5*(u[i,j,k]+u[i,j-1,k])
                
                UVy = (uavg2*vavg2 - uavg1*vavg1)/(ye[j]-ye[j-1])
                
                #UWz
                wavg2 = 0.5*(w[i,j,k]+w[i+1,j,k])
                wavg1 = 0.5*(w[i,j,k-1]+w[i+1,j,k-1])
                uavg2 = 0.5*(u[i,j,k]+u[i,j,k+1])
                uavg1 = 0.5*(u[i,j,k]+u[i,j,k-1])
                
                UWz = (uavg2*wavg2 - uavg1*wavg1)/(ze[k]-ze[k-1])
                
                #Ux2
                Ux2 = ((u[i+1,j,k]-u[i,j,k])/(xe[i+1]-xe[i]) - \
                       (u[i,j,k]-u[i-1,j,k])/(xe[i]-xe[i-1]))/ \
                       (xc[i+1]-xc[i])
                
                #Uy2
                Uy2 = ((u[i,j+1,k]-u[i,j,k])/(yc[j+1]-yc[j]) - \
                       (u[i,j,k]-u[i,j-1,k])/(yc[j]-yc[j-1]))/ \
                       (ye[j]-ye[j-1])
                
                Eu[i,j,k] = U2x + UVy + UWz - (Ux2 + Uy2)
    
    for i in range(1,Nx):
        for k in range(1,Nz):
            for j in range(1,Ny-1):
                # V2y
                vavg2 = 0.5*(v[i,j+1,k]+v[i,j,k])
                vavg1 = 0.5*(v[i,j,k]+v[i,j-1,k])
                
                V2y = ((vavg2)**2-(vavg1)**2)/(yc[j+1]-yc[j])
                
                #VUx
                uavg2 = 0.5*(u[i,j+1,k]+u[i,j,k])
                uavg1 = 0.5*(u[i-1,j,k]+u[i-1,j+1,k])
                vavg2 = 0.5*(v[i,j,k]+v[i+1,j,k])
                vavg1 = 0.5*(v[i,j,k]+v[i-1,j,k])
                
                VUx = (uavg2*vavg2 - uavg1*vavg1)/(xe[i]-xe[i-1])
                
                #VWz
                wavg2 = 0.5*(w[i,j,k]+w[i,j+1,k])
                wavg1 = 0.5*(w[i,j,k-1]+w[i,j+1,k-1])
                vavg2 = 0.5*(v[i,j,k]+v[i,j,k+1])
                vavg1 = 0.5*(v[i,j,k]+v[i,j,k-1])
                
                VWz = (vavg2*wavg2 - vavg1*wavg1)/(ze[k]-ze[k-1])
                
                #Vx2
                Vx2 = ((v[i+1,j,k]-v[i,j,k])/(xc[i+1]-xc[i]) - \
                       (v[i,j,k]-v[i-1,j,k])/(xc[i]-xc[i-1]))/ \
                       (xe[i]-xe[i-1])
                
                #Vy2
                Vy2 = ((v[i,j+1,k]-v[i,j,k])/(ye[j+1]-ye[j]) - \
                       (v[i,j,k]-v[i,j-1,k])/(ye[j]-ye[j-1]))/ \
                       (yc[j+1]-yc[j])
                
                Ev[i,j,k] = V2y + VUx + VWz - (Vx2 + Vy2)

For w-momentum equation, apart from the advection and the diffusion terms, the source term is also considered explicit.

    for j in range(1,Ny):
        for i in range(1,Nx):
            for k in range(1,Nz-1):
                # W2z
                wavg2 = 0.5*(w[i,j,k+1]+w[i,j,k])
                wavg1 = 0.5*(w[i,j,k]+w[i,j,k-1])
                
                W2z = ((wavg2)**2-(wavg1)**2)/(zc[k+1]-zc[k])
                
                #WVy
                vavg2 = 0.5*(v[i,j,k]+v[i,j,k+1])
                vavg1 = 0.5*(v[i,j-1,k]+v[i,j-1,k+1])
                wavg2 = 0.5*(w[i,j,k]+w[i,j+1,k])
                wavg1 = 0.5*(w[i,j,k]+w[i,j-1,k])
                
                WVy = (wavg2*vavg2 - wavg1*vavg1)/(ye[j]-ye[j-1])
                
                #WUx
                wavg2 = 0.5*(w[i+1,j,k]+w[i,j,k])
                wavg1 = 0.5*(w[i,j,k]+w[i-1,j,k])
                uavg2 = 0.5*(u[i,j,k]+u[i,j,k+1])
                uavg1 = 0.5*(u[i-1,j,k]+u[i-1,j,k+1])
                
                WUx = (uavg2*wavg2 - uavg1*wavg1)/(xe[i]-xe[i-1])
                
                #Wx2
                Wx2 = ((w[i+1,j,k]-w[i,j,k])/(xc[i+1]-xc[i]) - \
                       (w[i,j,k]-w[i-1,j,k])/(xc[i]-xc[i-1]))/ \
                       (xe[i]-xe[i-1])
                
                #Wy2
                Wy2 = ((w[i,j+1,k]-w[i,j,k])/(yc[j+1]-yc[j]) - \
                       (w[i,j,k]-w[i,j-1,k])/(yc[j]-yc[j-1]))/ \
                       (ye[j]-ye[j-1])
                
                Ew[i,j,k] = W2z + WVy + WUx - (Wx2 + Wy2) - (Ra/Pr)*T[i,j,k]

For the energy equation, the advection and diffusion terms are considered explicit.

    for i in range(1,Nx):
        for j in range(1,Ny):
            for k in range(1,Nz):
                #UTx
                Tavg2 = 0.5*(T[i,j,k]+T[i+1,j,k])
                Tavg1 = 0.5*(T[i-1,j,k]+T[i,j,k])
                
                UTx = (Tavg2*u[i,j,k]-Tavg1*u[i-1,j,k])/(xe[i]-xe[i-1])
                
                #VTy
                Tavg2 = 0.5*(T[i,j,k]+T[i,j+1,k])
                Tavg1 = 0.5*(T[i,j-1,k]+T[i,j,k])
                
                VTy = (Tavg2*v[i,j,k]-Tavg1*v[i,j-1,k])/(ye[j]-ye[j-1])
                
                #WTz
                Tavg2 = 0.5*(T[i,j,k]+T[i,j,k+1])
                Tavg1 = 0.5*(T[i,j,k-1]+T[i,j,k])
                
                WTz = (Tavg2*w[i,j,k]-Tavg1*w[i,j,k-1])/(ze[k]-ze[k-1])
                
                #Tx2
                Tx2 = ((T[i+1,j,k]-T[i,j,k])/(xc[i+1]-xc[i]) - \
                       (T[i,j,k]-T[i-1,j,k])/(xc[i]-xc[i-1]))/ \
                       (xe[i]-xe[i-1])
                
                #Ty2
                Ty2 = ((T[i,j+1,k]-T[i,j,k])/(yc[j+1]-yc[j]) - \
                       (T[i,j,k]-T[i,j-1,k])/(yc[j]-yc[j-1]))/ \
                       (ye[j]-ye[j-1])
                
                ET[i,j,k] = UTx + VTy + WTz - (1.0/Pr)*(Tx2 + Ty2)

## Calculation of the timestep

For the calculation of the timestep, CFL criterion and the diffusive limit are used.

    umax = np.max(u)
    vmax = np.max(v)
    wmax = np.max(w)
    delt = []
    c = 1.2
    for value, h in zip([umax,vmax,wmax],[hx,hy,hz]):
        if value != 0:
            delt.append(c*h/value)
    delt.append(0.5*c*hx**2)
    delt.append(0.5*c*hy**2)
    return np.min(delt)

## Runge-Kutta scheme formulation

In every RK substep, a velocity is predicted without considering the effect of pressure gradient. Due to the presence of implicit terms, a tri-diagonal matrix system is solved in order to obtain the predicted velocity values. Then, the predicted velocities are used to formulate the Pressure Poisson Equation. This equation is then solved to obtain the pressure value at every point of the pressure grid. Using these pressure values, the velocities are corrected, in order to satisfy the continuity equation. Advancement of temperature is simple. It is similar to the velocity prediction step. Temperature values are obtained by solving a tri-diagonal matrix system. Tri-diagonal matrix systems are solved using Thomas Algorithm.

    def rk(c1, c2, c3, delt, u, v, w, p, T, qu, qv, qw, qT, Eu, Ev, Ew, ET):
        deltrk = c3*delt
        qu = c1*qu + delt*Eu
        qv = c1*qv + delt*Ev
        qw = c1*qw + delt*Ew
        qT = c1*qT + delt*ET
        rhsu = u - c2*qu
        rhsv = v - c2*qv
        rhsw = w - c2*qw
        rhsT = T - c2*qT
        ustar = np.zeros((Nx,Ny+1,Nz+1))
        vstar = np.zeros((Nx+1,Ny,Nz+1))
        wstar = np.zeros((Nx+1,Ny+1,Nz))
        pstar = np.zeros((Nx+1,Ny+1,Nz+1))
        Tstar = np.zeros((Nx+1,Ny+1,Nz+1))
        ustar, vstar, wstar, pstar, Tstar = equatebc(ustar, \
                                        vstar, wstar, pstar, Tstar, \
                                        u, v, w, p, T)
        for i in range(1,Nx-1):
            for j in range(1,Ny):
                rtemp = rhsu[i,j,1:Nz].copy()
                rtemp[0] += deltrk/(hz**2)*ustar[i,j,0]
                rtemp[-1] += deltrk/(hz**2)*ustar[i,j,Nz]
                lower = -deltrk/(hz**2)*np.ones(Nz-2)
                upper = -deltrk/(hz**2)*np.ones(Nz-2)
                diag = (1+2*deltrk/(hz**2))*np.ones(Nz-1)
                utemp = ustar[i,j,1:Nz].copy()
                # Solving tri-diagonal matrix with Thomas algorithm
                utemp = thomas(lower, diag, upper, rtemp, utemp)
                ustar[i,j,1:Nz] = utemp.copy()
        for i in range(1,Nx):
            for j in range(1,Ny-1):
                rtemp = rhsv[i,j,1:Nz].copy()
                rtemp[0] += deltrk/(hz**2)*vstar[i,j,0]
                rtemp[-1] += deltrk/(hz**2)*vstar[i,j,-1]
                lower = -deltrk/(hz**2)*np.ones(Nz-2)
                upper = -deltrk/(hz**2)*np.ones(Nz-2)
                diag = (1+2*deltrk/(hz**2))*np.ones(Nz-1)
                vtemp = vstar[i,j,1:Nz].copy()
                # Solving tri-diagonal matrix with Thomas algorithm
                vtemp = thomas(lower, diag, upper, rtemp, vtemp)
                vstar[i,j,1:Nz] = vtemp.copy()
        for i in range(1,Nx):
            for j in range(1,Ny):
                rtemp = rhsT[i,j,1:Nz].copy()
                rtemp[0] += deltrk/(Pr*hz**2)*Tstar[i,j,0]
                rtemp[-1] += deltrk/(Pr*hz**2)*Tstar[i,j,-1]
                lower = -deltrk/(Pr*hz**2)*np.ones(Nz-2)
                upper = -deltrk/(Pr*hz**2)*np.ones(Nz-2)
                diag = (1+2*deltrk/(Pr*hz**2))*np.ones(Nz-1)
                Ttemp = Tstar[i,j,1:Nz].copy()
                # Solving tri-diagonal matrix with Thomas algorithm
                Ttemp = thomas(lower, diag, upper, rtemp, Ttemp)
                Tstar[i,j,1:Nz] = Ttemp.copy()
        for i in range(1,Nx):
            for j in range(1,Ny):
                rtemp = rhsw[i,j,1:Nz-1].copy()
                rtemp[0] += deltrk/(hz**2)*wstar[i,j,0]
                rtemp[-1] += deltrk/(hz**2)*wstar[i,j,-1]
                lower = -deltrk/(hz**2)*np.ones(Nz-3)
                upper = -deltrk/(hz**2)*np.ones(Nz-3)
                diag = (1+2*deltrk/(hz**2))*np.ones(Nz-2)
                wtemp = wstar[i,j,1:Nz-1].copy()
                # Solving tri-diagonal matrix with Thomas algorithm
                wtemp = thomas(lower, diag, upper, rtemp, wtemp)
                wstar[i,j,1:Nz-1] = wtemp.copy()
        div = np.zeros((Nx+1,Ny+1,Nz+1))
        div = divergence(ustar,vstar,wstar,div)
        rhsp = (1/deltrk)*div
        #Pressure Poisson Equation
        pstar = solve_p2(rhsp, pstar)
        #Correction of velocities
        px, py, pz = delp(pstar)
        ustar = ustar - deltrk*px
        vstar = vstar - deltrk*py
        wstar = wstar - deltrk*pz
        return ustar, vstar, wstar, pstar, Tstar, qu, qv, qw, qT

## Temporal Advancement

Each iteration involves three RK substeps, with the update of boundary conditions and calculation of explicit terms in between.

    while iteration < nt:
        dt = timestep(u,v,w)
        t += dt
        iteration += 1
        u_old = u.copy()
        v_old = v.copy()
        w_old = w.copy()
        T_old = T.copy()
        Eu, Ev, Ew, ET = explicit(u,v,w,T)
        qu = dt*Eu
        qv = dt*Ev
        qw = dt*Ew
        qT = dt*ET
        u, v, w, p, T, qu, qv, qw, qT = rk(0, 1/3, 1/3, dt, u, v, w, p, T, \
                                       qu, qv, qw, qT, Eu, Ev, Ew, ET)
        u, v, w, p, T = updateBC(u,v,w,p,T)
        Eu, Ev, Ew, ET = explicit(u,v,w,T)
        u, v, w, p, T, qu, qv, qw, qT = rk(-5/9, 15/16, 5/12, dt, u, v, w, p, T, \
                                       qu, qv, qw, qT, Eu, Ev, Ew, ET)
        u, v, w, p, T = updateBC(u,v,w,p,T)
        Eu, Ev, Ew, ET = explicit(u,v,w,T)
        u, v, w, p, T, qu, qv, qw, qT = rk(-153/128, 8/15, 1/4, dt, u, v, w, p, T, \
                                       qu, qv, qw, qT, Eu, Ev, Ew, ET)
        u, v, w, p, T = updateBC(u,v,w,p,T)
