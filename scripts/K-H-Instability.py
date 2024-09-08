"""
Python script which uses the Dedalus package to simulate continuously stratified shear instability
with a stable initial density profile, and consideration of viscous and diffusive effects.

Aspect ratio, resolution, Reynolds number, Prandtl number, and Richardson number are
all user specified arguments.

To run the script with more than 1 core and plot density:
    $ mpiexec -n <cores> python3 K-H-Instability.py...
    $ mpiexec -n <cores> python3 -m dedalus merge_procs snapshots --cleanup
    $ mpiexec -n <cores> python3 plot_slices.py snapshots/*.h5

Usage:
    K-H-Instability.py [options]

Options:
    --Re=<Reynolds>    Initial Reynolds number [default: 1965]
    --Pr=<Prandtl>     Prandtl number [default: 7]
    --Ri=<Richardson>  Richardson number [default: 0.08]
    --S=<scale>        Scale parameter- shear layer thickness to density layer thickness [default: 0.85]
    --aspect=<gamma>   Aspect ratio length:height [default: 2]
    --len=<length>     Horizontal length of the domain [default: 4]
    --res=<res>        Resolution in the veritcal. Horizontal resolution is given by aspect*res [default: 96]
    --dt=<dt>          Starting timestep size [default: 0.01]
    --sim_time=<time>  Dimensionless simulation time [default: 10]
    --basepath=<path>  Basepath for output files
"""

from docopt import docopt
from dedalus import public as de
from dedalus.extras import flow_tools
import numpy as np
import matplotlib.pyplot as plt
import h5py
from mpi4py import MPI
import time
import pathlib

import logging
logger = logging.getLogger(__name__)

args = docopt(__doc__)

len=float(args['--len'])
aspect=float(args['--aspect'])
res=int(args['--res'])

Lx, Lz = (len, len/aspect)
nx, nz = (int(aspect*res), res)

Pr = float(args['--Pr'])
Re = float(args['--Re'])
Ri = float(args['--Ri'])
S = float(args['--S'])

# construct a problem domain
x_basis = de.Fourier('x', nx, interval=(0, Lx), dealias=3/2)
z_basis = de.Chebyshev('z', nz, interval=(-Lz/2, Lz/2), dealias=3/2)
domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)

# set up IVP with BCs
problem = de.IVP(domain, variables=['u', 'w', 'p', 'rho', 'uz', 'wz', 'rhoz'])

problem.meta['w', 'uz', 'wz', 'rhoz']['z']['dirichlet'] = True

problem.parameters['Re'] = Re
problem.parameters['Pr'] = Pr
problem.parameters['Ri'] = Ri
problem.parameters['S'] = S

problem.parameters['kx'] = 1/Lx
problem.parameters['kz'] = 1/Lz

problem.add_equation("dt(u) + dx(p) - 1/Re*(dx(dx(u)) + dz(uz))          = -u*dx(u) - w*uz")
problem.add_equation("dt(w) + dz(p) - 1/Re*(dx(dx(u)) + dz(wz)) + Ri*rho = -u*dx(w) - w*wz")
problem.add_equation("dx(u) + wz                                         = 0")
problem.add_equation("dt(rho) - 1/(Pr*Re)*(dx(dx(rho)) + dz(rhoz))   = -u*dx(rho) - w*rhoz")

problem.add_equation("uz-dz(u)       = 0")
problem.add_equation("wz-dz(w)       = 0")
problem.add_equation("rhoz - dz(rho) = 0")

# free-slip boundaries
problem.add_bc("left(w)  = 0")
problem.add_bc("right(w) = 0", condition="(nx != 0)")
problem.add_bc("left(p) = 0",  condition="(nx == 0)")

problem.add_bc("left(uz)  = 0")
problem.add_bc("right(uz) = 0")

problem.add_bc("left(rhoz) = 0")
problem.add_bc("right(rhoz) = 0")

# build a solver
ts = de.timesteppers.RK443
solver =  problem.build_solver(ts)
logger.info('Solver built')

# construct initial conditions- restart from file if checkpoint exists
if not pathlib.Path('restart.h5').exists():
    x = domain.grid(0)
    z = domain.grid(1)

    u = solver.state['u']
    w = solver.state['w']
    rho = solver.state['rho']
    uz = solver.state['uz']
    wz = solver.state['wz']
    rhoz = solver.state['rhoz']

    gshape = domain.dist.grid_layout.global_shape(scales=1)
    slices = domain.dist.grid_layout.slices(scales=1)
    rand = np.random.RandomState(seed=42)
    noise = rand.standard_normal(gshape)[slices]

    # background + perturbations damped at walls
    zb, zt = z_basis.interval
    pert =  1e-3 * noise * (zt - z) * (z - zb)

    u['g'] = 0.5*np.tanh(2*z)
    w['g'] = 0.1*np.sin(4*np.pi*(x-Lx)/Lx)*np.exp(-(z/0.2)**2)
    rho['g'] = -(0.5*np.tanh(2*z*S)+pert)
    u.differentiate("z", out=uz)
    w.differentiate("z", out=wz)
    rho.differentiate("z", out=rhoz)

    dt = float(args['--dt'])
    stop_sim_time = float(args['--sim_time'])
    fh_mode = 'overwrite'

else:

    # Restart
    write, last_dt = solver.load_state('restart.h5', -1)

    # Timestepping and output
    dt = last_dt
    stop_sim_time = float(args['--sim_time'])
    fh_mode = 'append'

solver.stop_sim_time = stop_sim_time

#construct CFL condition for timestepping
CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=10, safety=0.8, max_dt=0.125, threshold=0.05)
CFL.add_velocities(('u', 'w'))

#filehandler to capture density profile
snapshots_file = str(args['--basepath']) + '_snapshots'
snapshots = solver.evaluator.add_file_handler(snapshots_file, sim_dt=0.1, max_writes=100, mode=fh_mode)
snapshots.add_task('rho', name='density')
snapshots.add_task('uz-dx(w)', name='vorticity')

analysis_file = str(args['--basepath']) + '_analysis'
analysis = solver.evaluator.add_file_handler(analysis_file, sim_dt=0.5, max_writes=100, mode=fh_mode)
analysis.add_task('Re', name='Re')
analysis.add_task('Pr', name='Pr')
analysis.add_task('Ri', name='Ri')
analysis.add_task('S', name='S')
analysis.add_task('interp(kx*integ_x(u), x=0)', name='h_avg_u')
analysis.add_task('interp(0.5*kz*kx*integ_z(integ_x(u**2+w**2)), x=0, z=0)', name='vol_avg_k')
analysis.add_task('interp(0.5*kz*integ_z((kx*integ_x(u))**2), x=0, z=0)', name='k_bar')
analysis.add_task('interp(Ri*kz*kx*integ_z(integ_x(rho*z)), x=0, z=0)', name='potential_E')

# filehandler to capture system state (for restarts)
state_file = str(args['--basepath'])+'_state'
state = solver.evaluator.add_file_handler(state_file, sim_dt=5, max_writes=100, mode=fh_mode)
state.add_system(solver.state)

#timestepping loop
try:
    logger.info('Starting loop')
    start_time = time.time()
    while solver.proceed:
        dt = CFL.compute_dt()
        dt = solver.step(dt)
        if (solver.iteration-1) % 10 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
    logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))
