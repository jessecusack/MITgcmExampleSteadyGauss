import numpy as np
import matplotlib.pyplot as plt
import shutil
import os


def lininc(n, Dx, dx0):
    a = (Dx-n*dx0)*2./n/(n+1)
    dx = dx0 + np.arange(1., n+1., 1.)*a
    return dx


Fr = 0.06
H = 2000.
h0 = 500.
om = 2.*np.pi/12.42/3600.
N0 = 5.2e-3
u0 = Fr*N0*h0


outdir = '../runs/RunFr{:1.3f}_tides'.format(Fr)

if not os.path.exists(outdir):
    os.mkdir(outdir)
if not os.path.exists(outdir+'/figs'):
    os.mkdir(outdir+'/figs')

shutil.copy('./gendata.py', outdir)

# These must match ../code/SIZE.h
ny = 1
nx = 4*20
nz = 25

# y direction:
dy = 1000
# x direction
xt = 410e3

nmid = 50
dx0 = 300.
nleft = int((nx-nmid)/2)
print(nleft)
nright = int((nx-nmid)/2)
dx = np.zeros(nx)
dxleft = np.flipud(lininc(nleft, 200.e3, dx0))
dxright = lininc(nright, 200.e3, dx0)
dx[0:nleft] = dxleft
dx[(nleft):(nleft+nmid)] = dx0
dx[(nleft+nmid):] = dxright
x = np.cumsum(dx)
x = x - x[int(np.floor(nx/2))]

with open(outdir+"/delXvar.bin", "wb") as f:
    dx.tofile(f)

# plot
if True:
    plt.figure()
    plt.plot(x/1000., dx)
    plt.xlim([-10, 10])
    plt.savefig(outdir+'/figs/dx.pdf')

# topo
sigma = 4000.  # m

topo = 1500*np.exp(-x*x/(sigma**2)) - 1500 + h0
# topo = h0*exp(-x*x/(3000**2))
print(topo.shape)
topo[topo < 0.] = 0.
topo = -H + topo
topo[topo < -H] = -H

# plot
if True:
    plt.figure()
    plt.plot(x/1.e3, topo)
    # plt.xlim([-20.,20.])
    plt.savefig(outdir+'/figs/topo.pdf')


with open(outdir+"/topo.bin", "wb") as f:
    topo.tofile(f)

# dz:
# dz is from the surface down (right?).  Its saved as positive.

dz = np.zeros(nz) + H/nz

with open(outdir+"/delZvar.bin", "wb") as f:
    dz.tofile(f)


# temperature profile...
g = 9.8
alpha = 2e-4
T0 = 28 + np.cumsum(N0**2/g/alpha*(-dz))

with open(outdir+"/TRef.bin", "wb") as f:
    T0.tofile(f)


# save T0 over whole domain
TT0 = np.tile(T0[:, np.newaxis], (1, nx))
with open(outdir+"/T0.bin", "wb") as f:
    TT0.tofile(f)


z = np.cumsum(dz)
# plot:
if True:
    plt.figure()
    plt.plot(T0, z)
    plt.savefig(outdir+'/figs/TO.pdf')

# Forcing for boundaries
dt = 3720.
time = np.arange(0, 12.*3720., dt)
print(time/3600./12.4)
om = 2*np.pi/12.40/3600
uw = u0*np.sin(om*time)
ue = u0*np.sin(om*time)

# plot:
if True:
    plt.figure()
    plt.plot(time/3600./12.4, ue, label='Ue')
    plt.plot(time/3600/12.4, uw, label='Uw')
    plt.legend()
    plt.xlabel('t/T')
    plt.ylabel('Vel')
    plt.title('{}'.format(time[-1]))
    plt.savefig(outdir+'/figs/Vels.pdf')

# try time,nz,ny...
uwn = np.tile(uw[:, np.newaxis, np.newaxis], (1, nz, ny))
uen = np.tile(ue[:, np.newaxis, np.newaxis], (1, nz, ny))

with open(outdir+"/Ue.bin", "wb") as f:
    uen.tofile(f)

with open(outdir+"/Uw.bin", "wb") as f:
    uwn.tofile(f)

t = np.tile(T0[np.newaxis, :, np.newaxis], (time.shape[0], 1, ny))

print(t.shape)
with open(outdir+"/Te.bin", "wb") as f:
    t.tofile(f)

with open(outdir+"/Tw.bin", "wb") as f:
    t.tofile(f)

# Copy some other files
shutil.copy('data', outdir+'/data')
shutil.copy('eedata', outdir)
shutil.copy('data.kl10', outdir)
shutil.copy('data.mnc', outdir)
shutil.copy('data.obcs', outdir)
shutil.copy('data.diagnostics', outdir)
shutil.copy('data.pkg', outdir+'/data.pkg')
# also store these.  They are small and helpful to document what we did

for nm in {'input', 'code', 'analysis'}:
    to_path = outdir+'/'+nm
    if os.path.exists(to_path):
        shutil.rmtree(to_path)
    shutil.copytree('../'+nm, outdir+'/'+nm)
