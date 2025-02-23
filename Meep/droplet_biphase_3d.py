## 3d lc biphase droplet 
# run as: mpirun -np 4 python3 droplet_biphase_3d.py
# python3 droplet_biphase_3d.py > output.txt &

## Note: grid initialization may take up to 30 min

## get libraries
import meep as mp                       # works for serial as well as parallel execution (depends on conda environment mp vs. pmp)
import numpy as np                      # advanced math functions
import matplotlib.pyplot as plt         # plot stuff
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable     # fit colorbar size to plot
import os                               # OS-independent directories
import argparse                         # read input parameters for script
import h5py                             # read/write h5 files
from scipy.interpolate import RegularGridInterpolator  # interpolation of sampled data
import shutil                           # copy files


### functions

## read input parameters
## call script as droplet_biphase_3d.py --theta ${theta} --phi ${phi}
def read_input_parameters():
    parser = argparse.ArgumentParser(description='Simulate electric field tranmission through emulsion droplets with meep for different tilt angles')
    parser.add_argument('--theta', default=0.0*np.pi/2.0, type=float, metavar='tiltTheta',help='theta tilt angle: 0 <= theta <= pi')
    parser.add_argument('--phi', default=0.0*np.pi/4.0, type=float, metavar='tiltPhi',help='phi tilt angle: 0 <= phi <= pi/4')
    # parser.add_argument('--outputName', default='test', type=float, metavar='name',help='name prefix of output file')
    args = parser.parse_args()
    tiltTheta = args.theta			# read theta
    tiltPhi = args.phi			    # read phi
    print(f"tilt angle theta: {tiltTheta:.3f}")
    print(f"tilt angle phi: {tiltPhi:.3f}")
    return [tiltTheta,tiltPhi]

# returns np.array
def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    if any(v): #if not all zeros then 
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        return np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    else:
        return np.eye(3) #cross of all zeros only occurs on identical directions

## convert meep vector3 to numpy vector
def mp2np(mpVector):
    return np.array(mpVector)





## set parameters and read input parameters
[tiltTheta,tiltPhi] = read_input_parameters()

trialIndex=16
# 8 - shell, interpolated field (slow, init > 75min), short
# 9 - hedgehog (faster) with issue at origin!, short
# 10 - hedgehog (faster) without issue at origin, short
# 11 - hedgehog meep (10% slower) without issue at origin, short
# 12 - shell, interpolated field (slow, init > 75min), long
# 13 - hedgehog (faster) without issue at origin, long, Iy
# 14 - shell, interpolated field (slow, > 75min), long, Iy
# 15 - shell, analytical approximation (slow, init > 10min), long, Iy (12329.2665 s for 50; 9459.0974 s for 35)
# 16 - shell, analytical approximation (slow, init > 10min), long, Iy, nice plot


## 10 procs: 1167.23s for grid init
## 6 procs: xxx s for grid init


nCPUs = mp.count_processors()
output_dir = os.path.join(".", f"output_ana_droplet_theta_{tiltTheta:.2f}_phi_{tiltPhi:.2f}_cpus_{nCPUs}_trial_{trialIndex}")
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)
shutil.copyfile(__file__, os.path.join(output_dir,__file__))

eps_averaging = False              # use subpixel averaging (True, default) or not (False)
force_complex_fields = True        # use real (False, default) or complex (True) fields in simulation?
use_symmetry = False               # use symmetry for speedup (True) or not (False)

# simulationTime = 100   # full simulation; for (16,16,24) too long
simulationTime = 35   # full simulation
# simulationTime = 70   # full simulation
# simulationTime = 1     # test simulation
# simulationTime = 5     # test simulation


cell = mp.Vector3(16, 16, 24)    # numbers are in microns
# resolution = 40                  # for interpolation spatial resolution: 40 pixels/micron
resolution = 30                  # spatial resolution: 30 pixels/micron


complexFactor = int(force_complex_fields) + 1       # False - x1 (real), True - x2 (real+imag)
floatSize=8         # double datatype
nFields=12        # 3d E,D,H (9) + diagonal eps (3)
requiredMemoryGB = complexFactor*floatSize*cell.x*cell.y*cell.z*1e-9*resolution**3
print(f"Required memory (GB): {requiredMemoryGB:.3f}")

nLC=1.6         # effective constant
nLCe=1.7074
nLCo=1.5343
epsLCe=nLCe*nLCe
epsLCo=nLCo*nLCo
nFC=1.282
nH2O=1.33

rDroplet=3.5                # in um
tiltDirector3d = mp.Vector3(np.cos(tiltPhi)*np.sin(tiltTheta), np.sin(tiltPhi)*np.sin(tiltTheta), np.cos(tiltTheta))

wvl = 0.550        # lambda = 550 nm; spatial scaling is by default: 1 micron
freq = 1.0/wvl     # frequency f, c=1 in meep units
bndWidth = 1.0     # width of absorbing boundary layer
freqInMedium = freq/nH2O
periodInMedium = 1.0/freqInMedium

pxWvl=resolution*(wvl/nLCe)                 # [res] = px/cell, [wvl/nLCe] = cell/lambda (in highest index medium)
print(f"px/wvl in highest index medium: {pxWvl:.2f}",)
if pxWvl < 8:
    print(f"Warning: Resolution might be too low. Px/wvl is {pxWvl:.2f} < 8")

sourceField = mp.Ex
observeField = mp.Ey

## load position and director vector data {x,y,z,director} from hdf5 file
file="./positionDirectorRegularGrid.hdf5"
hfile = h5py.File(file,'r')
xSamples = np.array(hfile.get('xVals'))
ySamples = np.array(hfile.get('yVals'))
zSamples = np.array(hfile.get('zVals'))
directorSamples = np.array(hfile.get('directorVals'))
rDroplet = 3.5
[xSamples,ySamples,zSamples] = [rDroplet*item for item in [xSamples,ySamples,zSamples]]
hfile.close()

## interpolate stuff
directorInterpolant = RegularGridInterpolator((xSamples, ySamples, zSamples), directorSamples)


## permittivity tensor of uniaxial LC aligned in cartesian coordinate system (director points in z-direction)
## in numpy to speed up later calculations
## points in e_z direction
epsilonLC = np.array([
    [epsLCo, 0, 0],
    [0, epsLCo, 0],
    [0, 0, epsLCe]
])

rotaMatrix = rotation_matrix_from_vectors(mp2np(tiltDirector3d),np.array([0,0,-1]))      # from tilted to upright; 
rotaMatrixInv = rotation_matrix_from_vectors(np.array([0,0,-1]),mp2np(tiltDirector3d))

mediumLC = mp.Medium(epsilon_diag=mp.Vector3(epsLCo, epsLCo, epsLCe))
mediumFC = mp.Medium(index=nFC)
mediumH2O = mp.Medium(index=nH2O)

def biphase_droplet(p):
    if p.dot(tiltDirector3d) >= 0.0:
        return mediumLC
    else:
        return mediumFC
biphase_droplet.do_averaging = eps_averaging    ## turn on subpixel averaging

bcBottom = np.array([0, 0, -1])
def bcHemisphere(x, y, z):
    eps = 0.000000001
    norm = np.sqrt(x**2 + y**2 + eps)
    return np.array([ 
        (x*z)/norm, 
        (y*z)/norm, 
        -(rDroplet**2 - z**2)/norm
    ])

def directorAnalyticalApproximation(p0):
    x,y,z = p0
    eps = 0.000000001
    if np.sqrt(x**2 + y**2) < rDroplet:
        nonNormDirector = bcBottom + z/np.sqrt(rDroplet**2 - (x**2 + y**2) + eps)*(bcHemisphere(x, y, z) - bcBottom)
        return nonNormDirector / np.linalg.norm(nonNormDirector)
    else:
        return np.array([0,0,0])

def biphase_droplet_lc(p):
    if p.norm() < rDroplet:
        if p.dot(tiltDirector3d) >= 0.0:                  # in meep
            p0 = rotaMatrix @ mp2np(p)                    # in numpy; position p0 in untilted configuration corresponding to tilted position p
            # director = directorInterpolant(p0)[0]         # in numpy; get director from sampled mma solution and boundary condition
            director = directorAnalyticalApproximation(p0)  # in numpy; get director from sampled mma solution and boundary condition
            tiltedDirector = rotaMatrixInv @ director     # in numpy; rotate director into tilted configuration
            rotaMatrixDirector = rotation_matrix_from_vectors([0,0,-1],tiltedDirector)   # in numpy; get rotation matrix that rotates director from natural state into tilted position
            rotatedEpsilon = rotaMatrixDirector @ epsilonLC @ rotaMatrixDirector.transpose()
            lc_epsilon_diag = mp.Vector3(rotatedEpsilon[0,0], rotatedEpsilon[1,1], rotatedEpsilon[2,2])
            lc_epsilon_offdiag = mp.Vector3(rotatedEpsilon[1,0], rotatedEpsilon[2,0], rotatedEpsilon[2,1])
            return mp.Medium(epsilon_diag=lc_epsilon_diag, epsilon_offdiag=lc_epsilon_offdiag)   # transformed permittivity tensor by director orientation
        else:
            return mediumFC
    else:
        return mediumH2O
biphase_droplet_lc.do_averaging = eps_averaging


def biphase_droplet_lc_hedgehog(p):
    if p.norm() < rDroplet:
        if p.dot(tiltDirector3d) >= 0.0:                  # in meep
            p0 = rotaMatrix @ mp2np(p)                    # in numpy; position p0 in untilted configuration corresponding to tilted position p
            if np.linalg.norm(p0) > 0.00001:
                director = p0/np.linalg.norm(p0)                     # in numpy; get director from sampled mma solution and boundary condition
            else:
                director = np.array([0,0,1])
            tiltedDirector = rotaMatrixInv @ director     # in numpy; rotate director into tilted configuration
            rotaMatrixDirector = rotation_matrix_from_vectors(np.array([0,0,-1]),tiltedDirector)   # in numpy; get rotation matrix that rotates director from natural state into tilted position
            rotatedEpsilon = rotaMatrixDirector @ epsilonLC @ rotaMatrixDirector.transpose()
            lc_epsilon_diag = mp.Vector3(rotatedEpsilon[0,0], rotatedEpsilon[1,1], rotatedEpsilon[2,2])
            lc_epsilon_offdiag = mp.Vector3(rotatedEpsilon[1,0], rotatedEpsilon[2,0], rotatedEpsilon[2,1])
            return mp.Medium(epsilon_diag=lc_epsilon_diag, epsilon_offdiag=lc_epsilon_offdiag)   # transformed permittivity tensor by director orientation
        else:
            return mediumFC
    else:
        return mediumH2O
biphase_droplet_lc_hedgehog.do_averaging = eps_averaging




## run simulation

print(f'Propagation of electromagnetic waves through tilted complex emulsion droplet')

geometry = [
    mp.Sphere(
        radius=rDroplet,
        center=mp.Vector3(),
        # material=biphase_droplet              # no director dependence
        material=biphase_droplet_lc
        # material=biphase_droplet_lc_hedgehog    # analytical test case
    )
]


sources = [
    mp.Source(
        mp.ContinuousSource(frequency=freq,is_integrated=True),     # is_integrated=True must be set for sources extending into PMLs
        center=mp.Vector3(0.0,0.0,5.0),
        size=mp.Vector3(cell.x,cell.y,0.0),
        component=sourceField
    )
]


pml_layers = [mp.PML(bndWidth,direction=mp.Z)]  # for reflecting or periodic boundary conditions
# pml_layers = [mp.PML(bndWidth)]


## exploit the mirror symmetry in structure+source: (tilt in e_x direction -> mirror symmetry along y-axis)
# TODO: rotate initial fields and readout fields to account for rotating polarizer
## exploit mirror symmetry (instead of others for debugging droplet)
if use_symmetry:
    symmetries = [mp.Mirror(mp.Y)]
else:
    symmetries = None


filename_prefix = f"theta_{tiltTheta:.3f}_phi_{tiltPhi:.3f}_"
offset = mp.Vector3(0,0,-5)
sim = mp.Simulation(
    cell_size=cell,
    boundary_layers=pml_layers,
    geometry=geometry,
    symmetries=symmetries,
    geometry_center=offset,
    default_material=mediumH2O,
    sources=sources,
    resolution=resolution,
    k_point = mp.Vector3(0,0,0),        # periodic boundary conditions
    filename_prefix=filename_prefix,
    progress_interval=10,               # print out progress information every 60 s (instead of 4s, to keep logfile clean)
    force_complex_fields=force_complex_fields,
    #Courant=0.5*0.55,
    ## subpixel averaging options
    eps_averaging=eps_averaging,
    # subpixel_tol=1e-4,           # default: 1e-4
    # subpixel_maxeval=1e4,        # default: 1e4
)
sim.use_output_directory(output_dir)

zMeasurePlane = -0.5*cell.z + 2*bndWidth + offset.z
# dtPlot = 1.0
dtPlot = 0.02
## custom stepper only takes sim as argument, no custom other arguments
output_xz_plane = mp.Volume(center=offset, size=mp.Vector3(cell.x, 0.0, cell.z))
output_xy_plane = mp.Volume(mp.Vector3(0,0,zMeasurePlane), size=mp.Vector3(cell.x - 2 * bndWidth, cell.y - 2 * bndWidth, 0))
   
def step_plotCrossSection(sim):         # custom stepper only takes sim as argument, no custom other arguments; use wrapper
    
    ## get data from xy-plane
    (x,y,_,_)=sim.get_array_metadata(vol=output_xy_plane)
    x1,x2 = [x[i] for i in (0, -1)]
    y1,y2 = [y[i] for i in (0, -1)]
    
    eFieldX = sim.get_array(component=mp.Ex, vol=output_xy_plane, cmplx=True)
    eFieldY = sim.get_array(component=mp.Ey, vol=output_xy_plane, cmplx=True)
    exAbs=np.abs(eFieldX)
    eyAbs=np.abs(eFieldY)
    ix=np.square(exAbs)
    iy=np.square(eyAbs)
    
    ## get data from xz-plane for colorbars (not needed in > v1.25)
    ex_xz_r = np.real(sim.get_array(component=mp.Ex, vol=output_xz_plane, cmplx=True).flatten())
    ey_xz_r = np.real(sim.get_array(component=mp.Ey, vol=output_xz_plane, cmplx=True).flatten())
    epsTrace3d = sim.get_epsilon()    # vol=output_xz_plane does not work :(
    print(epsTrace3d.shape)
    epsTrace = epsTrace3d[:,round(0.5*cell.y*resolution),:]
    
    ## get time
    t = round(sim.round_time()/dtPlot)
    
    ## plot trace of epsilon tensor
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(10, 8))
    ax=axs[0,0]
    plt.sca(ax)
    sim.plot2D(output_plane=output_xz_plane)
    plt.axhline(y=zMeasurePlane, xmin=bndWidth/cell.x, xmax=1-bndWidth/cell.x, color='b', linestyle='--')
    ax.set_title('medium')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", pad=0.1,size="5%")
    fig.add_axes(cax)
    mpl.colorbar.ColorbarBase(cax, cmap='gray', orientation="vertical", norm=mpl.colors.Normalize(vmin=np.min(epsTrace), vmax=np.max(epsTrace)))
    
    ## get+plot efield in xz-cross section
    ax=axs[0,1]
    plt.sca(ax)
    sim.plot2D(output_plane=output_xz_plane, fields=sourceField)        # , field_parameters={'colorbar':True}
    plt.axhline(y=zMeasurePlane, xmin=bndWidth/cell.x, xmax=1-bndWidth/cell.x, color='b', linestyle='--')
    ax.set_title('Re $E_x$')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", pad=0.1,size="5%")
    fig.add_axes(cax)
    mpl.colorbar.ColorbarBase(cax, cmap='RdBu', orientation="vertical", norm=mpl.colors.Normalize(vmin=np.min(ex_xz_r), vmax=np.max(ex_xz_r)))
    
    ax=axs[0,2]
    plt.sca(ax)
    sim.plot2D(output_plane=output_xz_plane, fields=observeField)       # , field_parameters={'colorbar':True}
    plt.axhline(y=zMeasurePlane, xmin=bndWidth/cell.x, xmax=1-bndWidth/cell.x, color='b', linestyle='--')
    ax.set_title('Re $E_y$')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", pad=0.1,size="5%")
    fig.add_axes(cax)
    mpl.colorbar.ColorbarBase(cax, cmap='RdBu', orientation="vertical", norm=mpl.colors.Normalize(vmin=np.min(ey_xz_r), vmax=np.max(ey_xz_r)))
    
    ax=axs[1,0]
    im = ax.imshow(ix+iy, cmap='gray', origin='lower', interpolation='none', extent=[x1,x2,y2,y1])
    im.set_clim(0,1)
    ax.set_title('$|E_x|^2+|E_y|^2$')
    ax.set_aspect(1)
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im,cax=cax)
    
    ax=axs[1,1]
    im = ax.imshow(ix, cmap='gray', origin='lower', interpolation='none', extent=[x1,x2,y2,y1])
    im.set_clim(0,0.5)
    ax.set_title('$|E_x|^2$')
    ax.set_aspect(1)
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im,cax=cax)
    
    ax=axs[1,2]
    im = ax.imshow(iy, cmap='gray', origin='lower', interpolation='none', extent=[x1,x2,y2,y1])
    im.set_clim(0,0.5)
    ax.set_title('$|E_y|^2$')
    ax.set_aspect(1)
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im,cax=cax)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, filename_prefix + f'cross_section_Ex_Ey_{t:04d}.png'))
    plt.close()


## execute simulation run and plot at regular intervals
tStartIntensityAveraging = simulationTime-1*periodInMedium
sim.run(
    ## at beginning (could also be before of run function)
    # mp.at_beginning(mp.output_epsilon(frequency=freq)),       ## save trace of epsilon
    mp.at_beginning(step_plotCrossSection),
    ## at every time step
    mp.at_every(dtPlot,step_plotCrossSection),      
    ## average output over last period
    # mp.after_time(tStartIntensityAveraging,mp.in_volume(output_xy_plane, mp.at_every(0.1*periodInMedium,mp.output_efield_x))),   # electric x
    # mp.after_time(tStartIntensityAveraging,mp.in_volume(output_xy_plane, mp.at_every(0.1*periodInMedium,mp.output_efield_y))),   # electric y
    # mp.after_time(tStartIntensityAveraging,mp.in_volume(output_xy_plane, mp.at_every(0.1*periodInMedium,mp.output_dpwr))),       # electric-field density
    ## 3d e-field for later paraview 3d visualization (TODO: fix name overlap); TODO: create paraview animation of focused light in 3d
    # mp.at_end(mp.output_efield_x),
    # mp.at_end(mp.output_efield_y),
    until=simulationTime
)


## save full output in different directory (and avoid name clash)
sim.use_output_directory(os.path.join(output_dir,"3d_fields"))
mp.output_efield_x(sim)
mp.output_efield_y(sim)
# later: run h5tovtk output.h5



## only run this on a single processor

#     ## full 3d image - seems buggy
#     plt.figure(dpi=100)
#     sim.plot3D(fields=mp.Ey)
#     if mp.am_master(): 
#       plt.savefig(os.path.join(output_dir, filename_prefix + 'full_Ey.png'))
#     plt.close()
    
    ## TODO: load h5 files from end
    
    ## TODO: average h5 files over one period