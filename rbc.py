import numpy as np
from numba import jit
import timeit
import scipy.sparse as sp

@jit(nopython=True)
def thomas(a,b,c,d,x):
    n = len(d)
    for i in range(1,n):
        w = a[i-1]/b[i-1]
        b[i] = b[i] - w*c[i-1]
        d[i] = d[i] - w*d[i-1]
    x[n-1] = d[n-1]/b[n-1]
    for i in range(n-2,-1,-1):
        x[i] = (d[i] - c[i]*x[i+1])/b[i]
    return x

def timestep(u,v,w):
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
    
@jit(nopython=True, fastmath=True)
def updateBC(u,v,w,p,T):
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
    # For velocities, no-slip at z-min and all-fixed at z-max
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
    return u, v, w, p, T

@jit(nopython=True, fastmath=True)
def equatebc(ustar,vstar,wstar,pstar,Tstar,u,v,w,p,T):
    # Periodic in x
    for j in range(Ny+1):
        for k in range(Nz+1):
            ustar[0,j,k] = u[Nx-2,j,k]
            ustar[Nx-1,j,k] = u[1,j,k]
            pstar[0,j,k] = p[Nx-1,j,k]
            pstar[Nx,j,k] = p[1,j,k]
            Tstar[0,j,k] = T[Nx-1,j,k]
            Tstar[Nx,j,k] = T[1,j,k]
            if j != Ny:
                vstar[0,j,k] = v[Nx-1,j,k]
                vstar[Nx,j,k] = v[1,j,k]
            if k != Nz:
                wstar[0,j,k] = w[Nx-1,j,k]
                wstar[Nx,j,k] = w[1,j,k]
    # Periodic in y
    for i in range(Nx+1):
        for k in range(Nz+1):
            vstar[i,0,k] = v[i,Ny-2,k]
            vstar[i,Ny-1,k] = v[i,1,k]
            pstar[i,0,k] = p[i,Ny-1,k]
            pstar[i,Ny,k] = p[i,1,k]
            Tstar[i,0,k] = T[i,Ny-1,k]
            Tstar[i,Ny,k] = T[i,1,k]
            if i != Nx:
                ustar[i,0,k] = u[i,Ny-1,k]
                ustar[i,Ny,k] = u[i,1,k]
            if k != Nz:
                wstar[i,0,k] = w[i,Ny-1,k]
                wstar[i,Ny,k] = w[i,1,k]
    # For velocities, no-slip at z-min and all-fixed at z-max
    # For temperature, all-fixed at z-min and zmax
    for i in range(Nx+1):
        for j in range(Ny+1):
            wstar[i,j,0] = 0
            wstar[i,j,Nz-1] = 0
            pstar[i,j,0] = p[i,j,1]
            pstar[i,j,Nz] = p[i,j,Nz-1]
            Tstar[i,j,0] = 2 - T[i,j,1] # T(zmin) = 1
            Tstar[i,j,Nz] = - T[i,j,Nz-1] # T(zmax) = 0
            if i != Nx:
                ustar[i,j,0] = - u[i,j,1]
                ustar[i,j,Nz] = - u[i,j,Nz-1]
            if j != Ny:
                vstar[i,j,0] = - v[i,j,1]
                vstar[i,j,Nz] = - v[i,j,Nz-1]
    return ustar, vstar, wstar, pstar, Tstar

def generate_A(ijklevels):
    numLevels = ijklevels.shape[1]
    A = []    
    for i in range(numLevels):
        Nx = ijklevels[0][i]
        Ny = ijklevels[1][i]
        Nz = ijklevels[2][i]
        nx = Nx-2
        ny = Ny-2
        nz = Nz-2
        hx = Lx/(Nx-1)
        hy = Ly/(Ny-1)
        hz = Lz/(Nz-1)
        diagonalsx = [(hz**2/hx**2)*np.ones(nx-1), (-2*hz**2/hx**2)*np.ones(nx), \
                      (hz**2/hx**2)*np.ones(nx-1)]
        diagonalsy = [(hz**2/hy**2)*np.ones(ny-1), (-2*hz**2/hy**2)*np.ones(ny), \
                      (hz**2/hy**2)*np.ones(ny-1)]
        diagonalsz = [np.ones(nz-1), -2*np.ones(nz), np.ones(nz-1)]
        Ax = sp.diags(diagonalsx,[-1,0,1])
        Ay = sp.diags(diagonalsy,[-1,0,1])
        Az = sp.diags(diagonalsz,[-1,0,1])
        A.append(sp.kronsum(sp.kronsum(Ax,Ay),Az,format='coo'))
    return A

def generateAfreqAdiag(A, ijklevels):
    numLevels = ijklevels.shape[1]
    Afreq = []
    Adiag = []
    for ilevel in range(numLevels):
        freq = np.bincount(A[ilevel].row)
        diag = np.where(A[ilevel].row == A[ilevel].col)[0]
        Afreq.append(freq)
        Adiag.append(diag)
    return Afreq, Adiag

@jit(nopython=True, fastmath=True)
def explicit(u,v,w,T):
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
    
    return Eu, Ev, Ew, ET

def generate_b(ijklevels, phi, f):
    Nx = ijklevels[0][0]
    Ny = ijklevels[1][0]
    Nz = ijklevels[2][0]
    hx = Lx/(Nx-1)
    hy = Ly/(Ny-1)
    hz = Lz/(Nz-1)
    nx = Nx-2
    ny = Ny-2
    nz = Nz-2    
    neq = nx*ny*nz
    g = np.empty(neq)
    b = np.empty(neq)
    for i in range(1,nx+1):
        for j in range(1,ny+1): 
            g[(i-1)+nx*(j-1):neq:nx*ny] = - ((hz**2/hx**2)*phi[i-1,j,1:nz+1] + \
              (hz**2/hx**2)*phi[i+1,j,1:nz+1] + (hz**2/hy**2)*phi[i,j+1,1:nz+1] + \
              (hz**2/hy**2)*phi[i,j-1,1:nz+1] \
              + phi[i,j,:nz] + phi[i,j,2:])
            b[(i-1)+nx*(j-1):neq:nx*ny] = hz**2*f[i,j,1:nz+1] \
            + g[(i-1)+nx*(j-1):neq:nx*ny]
    return b

def restrict(xfe, nfx, nfy, nfz, memory, notInMemory = True):
    ncx = int((nfx+1)/2) - 1
    ncy = int((nfy+1)/2) - 1
    ncz = int((nfz+1)/2) - 1
    for i, value in enumerate(memory[0]):
        if(value == (nfx,nfy,nfz)):
            R3D = memory[1][i]
            notInMemory = False
            break
    if notInMemory:
        block = np.array([1, 2, 1])
        Rx = sp.lil_matrix((ncx,nfx))
        Ry = sp.lil_matrix((ncy,nfy))
        Rz = sp.lil_matrix((ncz,nfz))
        for i in range(ncx):
            Rx[i, 2*i] = 0.25*block[0]
            Rx[i, 2*i+1] = 0.25*block[1]
            Rx[i, 2*i+2] = 0.25*block[2]
        for j in range(ncy):
            Ry[j, 2*j] = 0.25*block[0]
            Ry[j, 2*j+1] = 0.25*block[1]
            Ry[j, 2*j+2] = 0.25*block[2]
        for k in range(ncz):
            Rz[k, 2*k] = 0.25*block[0]
            Rz[k, 2*k+1] = 0.25*block[1]
            Rz[k, 2*k+2] = 0.25*block[2]
        R3D = sp.kron(Rz,sp.kron(Ry,Rx))
        memory[0].append((nfx,nfy,nfz))
        memory[1].append(R3D)
    xce = R3D.dot(xfe)
    
    return xce, ncx, ncy, ncz

@jit(nopython=True)
def gauss_seidel_update_optimized(Arow, Acol, Adata, b, x, freq, diag):
    k = len(freq)
    count = 0
    for i in range(k):
        temp = 0
        for j in range(count,count+freq[i]):
            if j == diag[i]:
                temp += b[i]
            else:
                temp += - x[Acol[j]]*Adata[j]
        x[i] = 1/Adata[diag[i]]*temp
        count += freq[i]
    return x

def gseidel_optimized(A, b, x, freq, diag, maxiters=10000, tol=1e-3):
    iteration = 0
#    n = int(len(b)**0.5) + 2
    error = tol + 1
    neq = len(x)
    factor = -0.5/(hz**2/hx**2+hz**2/hy**2+1)
    while iteration < maxiters and error > tol:
        x = gauss_seidel_update_optimized(A.row, A.col, A.data,b,x,freq,diag)
        r = factor*(b - A.dot(x))
#        r = computeResidual(A.row, A.col, A.data,b,x,freq)
        norm_r = np.linalg.norm(r)
        error = np.sqrt(norm_r**2/neq)
        iteration += 1
#        time = timeit.default_timer()
#        print('CPU Time = '+ str(time) + ', Iteration = '+ str(iteration) + ', |Residual| = '+ str(error) + '\n')
#        print('Nodes = '+ str(n) + ', Iteration = '+ str(iteration) + ', |Residual| = '+ str(error) + '\n')
    return x

def prolong(xce, ncx, ncy, ncz, memory, notInMemory = True):
    nfx = (ncx+1)*2 - 1
    nfy = (ncy+1)*2 - 1
    nfz = (ncz+1)*2 - 1
    for i, value in enumerate(memory[0]):
        if(value == (ncx,ncy,ncz)):
            I3D = memory[1][i]
            notInMemory = False
            break
    if notInMemory:
        block = np.array([1, 2, 1])
        Ix = sp.lil_matrix((nfx,ncx))
        Iy = sp.lil_matrix((nfy,ncy))
        Iz = sp.lil_matrix((nfz,ncz))
        for i in range(ncx):
            Ix[2*i, i] = 0.5*block[0]
            Ix[2*i+1, i] = 0.5*block[1]
            Ix[2*i+2, i] = 0.5*block[2]
        for j in range(ncy):
            Iy[2*j, j] = 0.5*block[0]
            Iy[2*j+1, j] = 0.5*block[1]
            Iy[2*j+2, j] = 0.5*block[2]
        for k in range(ncz):
            Iz[2*k, k] = 0.5*block[0]
            Iz[2*k+1, k] = 0.5*block[1]
            Iz[2*k+2, k] = 0.5*block[2]
        I3D = sp.kron(Iz,sp.kron(Iy,Ix))
        memory[0].append((ncx,ncy,ncz))
        memory[1].append(I3D)
    xfe = I3D.dot(xce)
    
    return xfe, nfx, nfy, nfz

def mg_update_optimized(b,x,level):
    #Pre-Smoothing
#    time = timeit.default_timer()
#    print('Pre-Smoothing started: ' + str(level))
    if level == numLevels-2:
        x = gseidel_optimized(A[level], b, x, Afreq[level], Adiag[level], maxiters=10)
    elif level == numLevels-3:
        x = gseidel_optimized(A[level], b, x, Afreq[level], Adiag[level], maxiters=5)
    elif level == numLevels-4:
        x = gseidel_optimized(A[level], b, x, Afreq[level], Adiag[level], maxiters=2)
#    elif level == numLevels-5:
#        x = gseidel(A[level], b, x, maxiters=2) 
#    elif level == numLevels-6:
#        x = gseidel(A[level], b, x, maxiters=1)
    else:
        x = gseidel_optimized(A[level], b, x, Afreq[level], Adiag[level], maxiters=1)
#    print('Pre-Smoothing finished: ' + str(level) + ' ' + str(timeit.default_timer() - time))
    #compute residual
    nfx = ijkLevels[0][level] - 2
    nfy = ijkLevels[1][level] - 2
    nfz = ijkLevels[2][level] - 2 
#    rf = computeResidual(A[level].row, A[level].col, A[level].data, b, x, Afreq[level])
    rf = b - A[level].dot(x)
    #Restriction
#    time = timeit.default_timer()
#    print('Restriction started: ' + str(level))
#    ncx = ijkLevels[0][level+1] - 2
#    ncy = ijkLevels[1][level+1]
#    ncz = ijkLevels[2][level+1] - 2
#    rc, ncx, ncy, ncz = restrict2(rf, intp(nfx), intp(nfy), intp(nfz))
#    rc = restrict(rf, nfx, nfy, nfz, restrict_memory)
    rc, ncx, ncy, ncz = restrict(rf, nfx, nfy, nfz, restrict_memory)
#    print('Restriction finished: ' + str(level) + ' ' + str(timeit.default_timer() - time))
    eps = np.zeros(len(rc))
    if level == numLevels-2:
        eps = gseidel_optimized(A[level+1], rc, eps, Afreq[level+1], Adiag[level+1])
    else:
        eps = mg_update_optimized(rc, eps, level+1)
#        if level == numLevels-2:
#            eps = gseidel(A[level+1], rc, eps, maxiters=100)
#        elif level == numLevels-3:
#            eps = gseidel(A[level+1], rc, eps, maxiters=60)
#        elif level == numLevels-4:
#            eps = gseidel(A[level+1], rc, eps, maxiters=40)
#        else:
#            eps = gseidel(A[level+1], rc, eps, maxiters=10)
#        rc = vcycle_update(eps, rc, level+1)
#    print('Prolongation started: ' + str(level))
#    epsf = prolong(eps, ncx, ncy, ncz, prolong_memory)
    epsf, nfx, nfy, nfz = prolong(eps, ncx, ncy, ncz, prolong_memory)
#    print('Prolongation finished: ' + str(level) + ' ' + str(timeit.default_timer() - time))
    x = x + epsf
#    rf = computeResidual(A[level].row, A[level].col, A[level].data, b, x, Afreq[level])
    rf = b - A[level].dot(x)
#    rc = restrict(rf, nfx, nfy, nfz, restrict_memory)
    rc, ncx, ncy, ncz = restrict(rf, nfx, nfy, nfz, restrict_memory)
    if level == numLevels-2:
        eps = gseidel_optimized(A[level+1], rc, eps, Afreq[level+1], Adiag[level+1])
    else:
        eps = mg_update_optimized(rc, eps, level+1)
#    epsf = prolong(eps, ncx, ncy, ncz, prolong_memory)
    epsf, nfx, nfy, nfz = prolong(eps, ncx, ncy, ncz, prolong_memory)
    x = x + epsf
    #Post-smoothing
#    print('Post-Smoothing started: ' + str(level))
    if level == numLevels-2:
        x = gseidel_optimized(A[level], b, x, Afreq[level], Adiag[level], maxiters=10)
    elif level == numLevels-3:
        x = gseidel_optimized(A[level], b, x, Afreq[level], Adiag[level], maxiters=5)
    elif level == numLevels-4:
        x = gseidel_optimized(A[level], b, x, Afreq[level], Adiag[level], maxiters=2)
#    elif level == numLevels-5:
#        x = gseidel(A[level], b, x, maxiters=2)
#    elif level == numLevels-6:
#        x = gseidel(A[level], b, x, maxiters=1)
    else:
        x = gseidel_optimized(A[level], b, x, Afreq[level], Adiag[level], maxiters=1)
#    print('Post-Smoothing finished: ' + str(level) + ' ' + str(timeit.default_timer() - time))
    
    return x

def multigrid_optimized(A,b,xe,maxiter=10000,tol = 1e-3):
    iteration = 0
    neq = len(xe)
    factor = -0.5/(hz**2/hx**2+hz**2/hy**2+1)
    r = factor*(b - A[0].dot(xe))
#    r = computeResidual(A[0].row, A[0].col, A[0].data, b[0], xe, Afreq[0])
    norm_r = np.linalg.norm(r)
    error = np.sqrt(norm_r**2/neq)
#    time = timeit.default_timer()
    #f.write(str(time) + '\t' + str(iteration) + '\t' + str(norm_r) + '\n')
#    print('CPU Time = '+ str(time) + ', V-cycle Iteration = ' + str(iteration) + ', |Residual| = ' + str(error))

    while iteration < maxiter and error > tol:
#        print(timeit.default_timer())
        xe = mg_update_optimized(b, xe, 0)
#        print(timeit.default_timer())
        r = factor*(b - A[0].dot(xe))
#        r = computeResidual(A[0].row, A[0].col, A[0].data, b[0], xe, Afreq[0])    
        norm_r = np.linalg.norm(r)
        error = np.sqrt(norm_r**2/neq)
        iteration += 1
#        time = timeit.default_timer()
        #f.write(str(time) + '\t' + str(iteration) + '\t' + str(norm_r) + '\n')
#        print('CPU Time = '+ str(time) + ', V-cycle Iteration = ' + str(iteration) + ', |Residual| = ' + str(error))
    return xe

@jit(nopython=True)
def initguess(d,x):
    N = len(x)
    a = 0.0625*np.ones(N-1)
    b = -2.25*np.ones(N)
    c = 0.0625*np.ones(N-1)
    return thomas(a,b,c,d,x)

def solve_p2(rhs,p):
    b = generate_b(ijkLevels, p, rhs) 
    nx = Nx-1
    ny = Ny-1
    nz = Nz-1
    neq = nx*ny*nz
#    p1d = np.zeros(neq)
#    p1d= np.ones(neq)
#    p1d = -(1/6)*b[0].copy()
    p1d = np.empty(neq)
    p1d = initguess(b.copy(),p1d)
    p1d = multigrid_optimized(A,b,p1d)
    p[1:nx+1,1:ny+1,1:nz+1] = p1d.reshape((nx,ny,nz),order='F')
    return p

@jit(nopython=True, fastmath=True)
def divergence(u,v,w,div):
    for k in range(1,Nz):
        for j in range(1,Ny):
            for i in range(1,Nx):    
                div[i,j,k] = (u[i,j,k]-u[i-1,j,k])/(xe[i]-xe[i-1]) + \
                             (v[i,j,k]-v[i,j-1,k])/(ye[j]-ye[j-1]) + \
                             (w[i,j,k]-w[i,j,k-1])/(ze[k]-ze[k-1])
    
    return div

@jit(nopython=True, fastmath=True)
def delp(p):
    Px = np.zeros((Nx,Ny+1,Nz+1))
    Py = np.zeros((Nx+1,Ny,Nz+1))
    Pz = np.zeros((Nx+1,Ny+1,Nz))
    for k in range(1,Nz):
        for j in range(1,Ny):
            for i in range(1,Nx-1):
                Px[i,j,k] = (p[i+1,j,k] - p[i,j,k])/(xc[i+1] - xc[i])
    for i in range(1,Nx):
        for k in range(1,Nz):
            for j in range(1,Ny-1):
                Py[i,j,k] = (p[i,j+1,k] - p[i,j,k])/(yc[j+1] - yc[j])
    for j in range(1,Ny):
        for i in range(1,Nx):
            for k in range(1,Nz-1):
                Pz[i,j,k] = (p[i,j,k+1] - p[i,j,k])/(zc[k+1] - zc[k])
    return Px, Py, Pz

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
            wtemp = thomas(lower, diag, upper, rtemp, wtemp)
            wstar[i,j,1:Nz-1] = wtemp.copy()
#    ustar, vstar, wstar, pstar = updateBC(ustar, vstar, wstar, pstar)
    div = np.zeros((Nx+1,Ny+1,Nz+1))
    div = divergence(ustar,vstar,wstar,div)
    rhsp = (1/deltrk)*div
    pstar = solve_p2(rhsp, pstar)
    px, py, pz = delp(pstar)
    ustar = ustar - deltrk*px
    vstar = vstar - deltrk*py
    wstar = wstar - deltrk*pz
    return ustar, vstar, wstar, pstar, Tstar, qu, qv, qw, qT

def residual(Eu, Ev, Ew, ET, u, v, w, p, T, u_old, v_old, w_old, T_old, delt):
    Iu = np.zeros((Nx,Ny+1,Nz+1))
    Iv = np.zeros((Nx+1,Ny,Nz+1))
    Iw = np.zeros((Nx+1,Ny+1,Nz))
    IT = np.zeros((Nx+1,Ny+1,Nz+1))
    div = np.zeros((Nx+1,Ny+1,Nz+1))
    
    for k in range(1,Nz):
        for j in range(1,Ny):
            for i in range(1,Nx-1):
                Iu[i,j,k] = ((u[i,j,k+1]-u[i,j,k])/(zc[k+1]-zc[k]) - \
                            (u[i,j,k]-u[i,j,k-1])/(zc[k]-zc[k-1]))/ \
                            (ze[k]-ze[k-1])
    
    for i in range(1,Nx):
        for k in range(1,Nz):
            for j in range(1,Ny-1):
                Iv[i,j,k] =((v[i,j,k+1]-v[i,j,k])/(zc[k+1]-zc[k]) - \
                           (v[i,j,k]-v[i,j,k-1])/(zc[k]-zc[k-1]))/ \
                           (ze[k]-ze[k-1])
    
    for j in range(1,Ny):
        for i in range(1,Nx):
            for k in range(1,Nz-1):
                Iw[i,j,k] = ((w[i,j,k+1]-w[i,j,k])/(ze[k+1]-ze[k]) - \
                       (w[i,j,k]-w[i,j,k-1])/(ze[k]-ze[k-1]))/ \
                       (zc[k+1]-zc[k])
                       
    for i in range(1,Nx):
        for j in range(1,Ny):
            for k in range(1,Nz):
                IT[i,j,k] = (1/Pr)*((T[i,j,k+1]-T[i,j,k])/(zc[k+1]-zc[k]) - \
                       (T[i,j,k]-T[i,j,k-1])/(zc[k]-zc[k-1]))/ \
                       (ze[k]-ze[k-1])
    
    px, py, pz = delp(p)
    
    div = divergence(u,v,w,div)
    
    residualU = (u-u_old)/delt + Eu + px - Iu
    residualV = (v-v_old)/delt + Ev + py - Iv
    residualW = (w-w_old)/delt + Ew + pz - Iw
    residualT = (T-T_old)/delt + ET - IT
    
    normU = np.linalg.norm(residualU)
    normV = np.linalg.norm(residualV)
    normW = np.linalg.norm(residualW)
    normT = np.linalg.norm(residualT)
    normdiv = np.linalg.norm(div)
    
#    print(np.linalg.norm(Eu),np.linalg.norm(px),np.linalg.norm(Iu))
#    print(np.linalg.norm(Ev),np.linalg.norm(py),np.linalg.norm(Iv))
#    print(np.linalg.norm(Ew),np.linalg.norm(pz),np.linalg.norm(Iw))
#    print(normdiv)
    
    res = np.sqrt(1/((Nx)*(Ny+1)*(Nz+1))*normU**2+\
                  1/((Nx+1)*(Ny)*(Nz+1))*normV**2+\
                  1/((Nx+1)*(Ny+1)*(Nz))*normW**2+\
                  1/((Nx+1)*(Ny+1)*(Nz+1))*normT**2+\
                  1/((Nx+1)*(Ny+1)*(Nz+1))*normdiv**2)
    
    return res

@jit(nopython=True)
def convertToNodalData(u,v,w,p,T):
    Nx = u.shape[0]
    Ny = v.shape[1]
    Nz = w.shape[2]
    
    un = np.empty((Nx,Ny,Nz))
    vn = np.empty((Nx,Ny,Nz))
    wn = np.empty((Nx,Ny,Nz))
    pn = np.empty((Nx,Ny,Nz))
    Tn = np.empty((Nx,Ny,Nz))
    
    utemp = np.empty((Nx,Ny,Nz+1))
    vtemp = np.empty((Nx+1,Ny,Nz))
    wtemp = np.empty((Nx,Ny+1,Nz))
    
    for k in range(Nz+1):
        for i in range(Nx):
            utemp[i,:,k] = 0.5*(u[i,:-1,k]+u[i,1:,k])
    
    for j in range(Ny):
        for i in range(Nx):
            un[i,j,:] = 0.5*(utemp[i,j,:-1]+utemp[i,j,1:])

    for i in range(Nx+1):
        for j in range(Ny):
            vtemp[i,j,:] = 0.5*(v[i,j,:-1]+v[i,j,1:])
    
    for j in range(Ny):
        for k in range(Nz):
            vn[:,j,k] = 0.5*(vtemp[1:,j,k]+vtemp[:-1,j,k])

    for j in range(Ny+1):
        for k in range(Nz):
            wtemp[:,j,k] = 0.5*(w[:-1,j,k]+w[1:,j,k])
    
    for k in range(Nz):
        for i in range(Nx):
            wn[i,:,k] = 0.5*(wtemp[i,1:,k]+wtemp[i,:-1,k])
    
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                pn[i,j,k] = (1/8)*(p[i,j,k] + p[i+1,j,k] +\
                 p[i,j+1,k] + p[i,j,k+1] + p[i+1,j+1,k] + p[i,j+1,k+1] +\
                 p[i+1,j,k+1] + p[i+1,j+1,k+1])
    
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                Tn[i,j,k] = (1/8)*(T[i,j,k] + T[i+1,j,k] +\
                 T[i,j+1,k] + T[i,j,k+1] + T[i+1,j+1,k] + T[i,j+1,k+1] +\
                 T[i+1,j,k+1] + T[i+1,j+1,k+1])
            
    return un, vn, wn, pn, Tn

def savedata(x,y,z,u,v,w,p,T,imax,jmax,kmax,iters):
    file = open(r"C:\Users\debar\Desktop\python_tecplot\rbc65\rbc65_"+str(iters)+".dat",'w')
    file.write("TITLE = \"Rayleigh-Benard Convection Data\" \n")
    file.write("VARIABLES = \"X\" \"Y\" \"Z\" \"U\" \"V\" \"W\" \"P\" \"T\" \n")
    file.write("ZONE \n")
    file.write("ZONETYPE = ORDERED, I = " + str(imax) + \
               ", J = " + str(jmax) + ", K = " + str(kmax) + "\n")
    file.write("DATAPACKING = POINT \n")
    for k in range(kmax):
        for j in range(jmax):
            for i in range(imax):
                file.write(str(x[i,j,k]) + '\t' + str(y[i,j,k]) + '\t' + \
                           str(z[i,j,k]) + '\t' + str(u[i,j,k]) + '\t' + \
                           str(v[i,j,k]) + '\t' + str(w[i,j,k]) + '\t' + \
                           str(p[i,j,k]) + '\t' + str(T[i,j,k]) + '\n')
    file.close()

def saveprobedata(w,T,iters,time):
    probeW = open(r"C:\Users\debar\Desktop\python_tecplot\rbc65\rbc65_wdata.dat",'a')
    probeT = open(r"C:\Users\debar\Desktop\python_tecplot\rbc65\rbc65_Tdata.dat",'a')
    if iters == 0:
        probeW.write("W-Velocity data at position x = 4, y = 4, z = 0.5 \n")
        probeW.write("Iteration \t Time \t w-velocity \n")
        probeT.write("Temperature data at position x = 4, y = 4, z = 0.5 \n")
        probeT.write("Iteration \t Time \t Temperature \n")
    wprobe = 0.5*(w[int((Nx+1)/2),int((Ny+1)/2),int(Nz/2)] + \
                    w[int((Nx+1)/2),int((Ny+1)/2),int(Nz/2)+1])
    Tprobe = T[int((Nx+1)/2),int((Ny+1)/2),int((Nz+1)/2)]
    probeW.write(str(iters) + "\t" + str(time) + "\t" + str(wprobe) + "\n")
    probeT.write(str(iters) + "\t" + str(time) + "\t" + str(Tprobe) + "\n")
    probeW.close()
    probeT.close()
    
def saverestartdata(u,v,w,p,T,iters,time):
    file = open(r"C:\Users\debar\Desktop\python_tecplot\rbc65\rbc65res_"+str(iters)+".dat",'w')
    file.write(str(iters) + '\n')
    file.write(str(time) + '\n')
    
    for k in range(Nz+1):
        for j in range(Ny+1):
            for i in range(Nx):
                file.write(str(u[i,j,k])+ '\n')
    
    for k in range(Nz+1):
        for j in range(Ny):
            for i in range(Nx+1):
                file.write(str(v[i,j,k])+ '\n')
    
    for k in range(Nz):
        for j in range(Ny+1):
            for i in range(Nx+1):
                file.write(str(w[i,j,k])+ '\n')
    
    for k in range(Nz+1):
        for j in range(Ny+1):
            for i in range(Nx+1):
                file.write(str(p[i,j,k])+ '\n')
    
    for k in range(Nz+1):
        for j in range(Ny+1):
            for i in range(Nx+1):
                file.write(str(T[i,j,k])+ '\n')

def readdata(file):
    iters = int(file.readline())
    time = float(file.readline())
    a = np.loadtxt(file, dtype='float64')
    uend = Nx*(Ny+1)*(Nz+1)
    vend = uend + Ny*(Nx+1)*(Nz+1)
    wend = vend + Nz*(Nx+1)*(Ny+1)
    pend = wend + (Nx+1)*(Ny+1)*(Nz+1)
    u = a[:uend].reshape((Nx,Ny+1,Nz+1),order='F')
    v = a[uend:vend].reshape((Nx+1,Ny,Nz+1),order='F')
    w = a[vend:wend].reshape((Nx+1,Ny+1,Nz),order='F')
    p = a[wend:pend].reshape((Nx+1,Ny+1,Nz+1),order='F')
    T = a[pend:].reshape((Nx+1,Ny+1,Nz+1),order='F')
    return u, v, w, p, T, iters, time

time = timeit.default_timer()
#print(timeit.default_timer())    
Ra = 1e4
Pr = 1
nt = 1000
ntrest = 100
ntdsk = 50
ntprobe = 25
restart = False
Nx = 64
Ny = 64
Nz = 32
N = Nx*Ny*Nz
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
#Initialise variables
u = np.zeros((Nx,Ny+1,Nz+1))
v =  np.zeros((Nx+1,Ny,Nz+1))
w = np.zeros((Nx+1,Ny+1,Nz))
p = np.zeros((Nx+1,Ny+1,Nz+1))
T = np.zeros((Nx+1,Ny+1,Nz+1))
#boundary conditions
#print(timeit.default_timer())
u, v, w, p, T = updateBC(u,v,w,p,T)
t = 0
iteration = 0
#restart
if restart:
    re_file = open(r"C:\Users\debar\Desktop\python_tecplot\rbc65\rbc65res_10.dat",'r')
    u, v, w, p, T, iteration, t = readdata(re_file)
    re_file.close()
#probe
saveprobedata(w,T,iteration,t)

prolong_memory = [[],[]]
restrict_memory = [[],[]]
numLevels = 5
ijkLevels = np.array([[65, 33, 17, 9, 5],\
                      [65, 33, 17, 9, 5],\
                      [33, 17, 9, 5, 3]])
A = generate_A(ijkLevels)
Afreq, Adiag = generateAfreqAdiag(A, ijkLevels)

while iteration < nt:
#    print(timeit.default_timer())
    dt = timestep(u,v,w)
#    dt_2 = timestep_updated(u,v,w)
#    dt = 0.01
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
    error = residual(Eu, Ev, Ew, ET, u, v, w, p, T, u_old, v_old, w_old, T_old, dt)
#    res_file = open("/scratch/debartha/advcfd/ldc/grid3/residual.txt",'a')
#    res_file.write('Iteration = ' + str(iteration) + ', time = ' + str(t) + \
#          ', CPU time = ' + str(timeit.default_timer()-time) + ', |Residual| = '\
#          + str(error) + '\n')
#    res_file.close()
    print('Iteration = ' + str(iteration) + ', time = ' + str(t) + \
          ', CPU time = ' + str(timeit.default_timer()-time) + ', |Residual| = '\
          + str(error))
#    print('dt = ' + str(dt) + ' , dt2 = ' + str(dt_2))
    if iteration % ntdsk == 0:
        un, vn, wn, pn, Tn = convertToNodalData(u,v,w,p,T)
        savedata(x,y,z,un,vn,wn,pn,Tn,Nx,Ny,Nz,iteration)
    if iteration % ntrest == 0:
        saverestartdata(u,v,w,p,T,iteration,t)
    if iteration % ntprobe == 0:
        saveprobedata(w,T,iteration,t)