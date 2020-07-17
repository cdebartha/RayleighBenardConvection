# Rayleigh Benard Convection

## Governing Equations

Incompressible Navier-Stokes equations with body force <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/930600685a1f78702b6e7cc965660e2d.svg?invert_in_darkmode" align=middle width=35.799139199999985pt height=32.42016360000002pt/> :

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/719f23136027ff8435fc8451db709992.svg?invert_in_darkmode" align=middle width=58.90396379999999pt height=11.232861749999998pt/></p> (1)

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/e2524fbcd54b0d0bc01df9da60a4d79a.svg?invert_in_darkmode" align=middle width=264.78570525pt height=33.81208709999999pt/></p> (2)

where, <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/3d96202da71b051caff14815ac3573bd.svg?invert_in_darkmode" align=middle width=89.99601764999998pt height=24.65753399999998pt/> are the velocities in <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/86ad5e70a1497f100df341d42edabbdd.svg?invert_in_darkmode" align=middle width=30.82954049999999pt height=14.15524440000002pt/> and <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/f93ce33e511096ed626b4719d50f17d2.svg?invert_in_darkmode" align=middle width=8.367621899999993pt height=14.15524440000002pt/> directions respectively, <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/2ec6e630f199f589a2402fdf3e0289d5.svg?invert_in_darkmode" align=middle width=8.270567249999992pt height=14.15524440000002pt/> is the pressure, <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?invert_in_darkmode" align=middle width=8.49888434999999pt height=14.15524440000002pt/> is the density and <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/07617f9d8fe48b4a7b3f523d6730eef0.svg?invert_in_darkmode" align=middle width=9.90492359999999pt height=14.15524440000002pt/> is the dynamic viscosity.

## Boussinesq Approximation

This approximation is based on the fact that variation of density induced by the variation of temperature can be everywhere except in the buoyancy term.

Say, the density is expressed as,

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/dcd7ec0bd28b17c9fa06026d90e628b6.svg?invert_in_darkmode" align=middle width=163.83224879999997pt height=16.438356pt/></p> (3)

where, <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/4f649235e80d795f0a153fc88d1e7d44.svg?invert_in_darkmode" align=middle width=14.373536099999988pt height=14.15524440000002pt/> and <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/f3087bed09ac6d957af8ecfe91212479.svg?invert_in_darkmode" align=middle width=15.480837899999988pt height=22.465723500000017pt/> are the density and the temperature of the fluid at the cold plate respectively.

Substituting <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?invert_in_darkmode" align=middle width=8.49888434999999pt height=14.15524440000002pt/> in the momentum equation, we get,

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/de1fee5091fd8a35af2e25cf18e6618b.svg?invert_in_darkmode" align=middle width=535.8008238pt height=39.452455349999994pt/></p> (4)

Now, applying the Boussinesq approximation, we neglect the dependence of density in the inertia terms.

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/2fc49f8455f21ed28b23b830766e1fcb.svg?invert_in_darkmode" align=middle width=417.58055955pt height=39.452455349999994pt/></p>
<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/62a04feb41e0b6fc9b7dca0ac4626f97.svg?invert_in_darkmode" align=middle width=452.56479015pt height=39.452455349999994pt/></p> (5)

Define:

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/5f4c1b3e7092df86604873fe590e1e43.svg?invert_in_darkmode" align=middle width=95.81168355pt height=19.5105966pt/></p> (6)

Therefore,

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/8d922d8b36f7361fa3e46966a3b41efb.svg?invert_in_darkmode" align=middle width=129.4247526pt height=19.5105966pt/></p> (7)

Substituting (7) in (5), we get,

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/16dad1c6b6914e6caa2bc032d29281ec.svg?invert_in_darkmode" align=middle width=385.40603805pt height=39.452455349999994pt/></p> (8)

Since, the number of variables is more than the number of
equations, we need energy equation to find <img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/2f118ee06d05f3c2d98361d9c30e38ce.svg?invert_in_darkmode" align=middle width=11.889314249999991pt height=22.465723500000017pt/>.

<p align="center"><img src="https://raw.githubusercontent.com/cdebartha/RayleighBenardConvection/master/svgs/3a816856923535f3eae73abaeda2e21e.svg?invert_in_darkmode" align=middle width=148.6011054pt height=33.81208709999999pt/></p> (9)

