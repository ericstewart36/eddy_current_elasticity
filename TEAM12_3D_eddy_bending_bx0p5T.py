# %% [markdown]
# #  TEAM problem 12, 3D
# 
# - Code for large deformation magneto-elasticity of conducting materials. 
# 
# -  We use 3-D magnetic vector potential $\mathbf{a}_{\text{\tiny R}}$ as the magnetic degree of freedom, with 
# $\mathbf{b}_\text{\tiny R} = \text{Curl } \mathbf{a}_\text{\tiny R}$.
# 

# %% [markdown]
# ### Units
# 
# - Basic:
#     - Length: mm
#     -   Time: s
#     -   Mass: kg
#     - Charge: C
#  
# - Derived: 
#     - Force: mN
#     - Pressure: kPa
#     - Energy: $\mu\text{J}$
#     - Current: A
#     - Mag. flux density: T
#     - Electrical conductivity: MS/mm, such that 1 MS/m = $10^{-3}$ MS/mm   
# 
# ### Software:
# - Dolfinx v0.9.0
# 
# Eric M. Stewart
# 
# stewaei@ucmail.uc.edu
# 
# Summer 2025
# 

# %% [markdown]
# # Import modules

# %%
# Import FEnicSx/dolfinx
import dolfinx

# For numerical arrays
import numpy as np

# For MPI-based parallelization
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# PETSc solvers
from petsc4py import PETSc

# specific functions from dolfinx modules
from dolfinx import fem, mesh, io, plot, log
from dolfinx.fem import (Constant, dirichletbc, Function, functionspace, Expression)
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.io import VTXWriter, XDMFFile


# specific functions from ufl modules
import ufl
from ufl import (TestFunctions, TestFunction, TrialFunction, Identity, grad, det, div, dev, inv, tr, sqrt, conditional ,\
                 gt, dx, inner, derivative, dot, ln, split, exp, eq, cos, sin, acos, ge, le, outer, tanh, cosh, atan, atan2)

# basix finite elements
import basix
from basix.ufl import element, mixed_element, quadrature_element

# Matplotlib for plotting
import matplotlib.pyplot as plt
plt.close('all')

# For timing the code
from datetime import datetime

# this forces the program to still print (but only from one CPU) 
# when run in parallel.
def mprint(*argv):
    if rank==0:
        out=""
        for arg in argv:
            out = out+ str(argv)
        print(out, flush=True)

# Set level of detail for log messages (integer)
# Guide:
# CRITICAL  = 50, // errors that may lead to data corruption
# ERROR     = 40, // things that HAVE gone wrong
# WARNING   = 30, // things that MAY go wrong later
# INFO      = 20, // information of general interest (includes solver info)
# PROGRESS  = 16, // what's happening (broadly)
# TRACE     = 13, // what's happening (in detail)
# DBG       = 10  // sundry
#
log.set_log_level(log.LogLevel.WARNING)

# %% [markdown]
# # Define geometry

# %%

# Read in the 3D mesh and cell tags
with XDMFFile(MPI.COMM_WORLD,"meshes/TEAM_12_v3.xdmf",'r') as infile:
    domain = infile.read_mesh(name="Grid",xpath="/Xdmf/Domain")
    cell_tags = infile.read_meshtags(domain,name="Grid")
domain.topology.create_connectivity(domain.topology.dim, domain.topology.dim-1)

# Also read in 2D facets for applying BCs
with XDMFFile(MPI.COMM_WORLD,"meshes/facet_TEAM_12_v3.xdmf",'r') as infile:
    facet_tags = infile.read_meshtags(domain,name="Grid")

x = ufl.SpatialCoordinate(domain)

# %%
# Define the boundary integration measure "ds" using the facet tags,
# also specify the number of surface quadrature points.
ds = ufl.Measure('dS', domain=domain, subdomain_data=facet_tags, metadata={'quadrature_degree':2})

# Define the volume integration measure "dx" 
# also specify the number of volume quadrature points.
dx = ufl.Measure('dx', domain=domain, subdomain_data=cell_tags, metadata={'quadrature_degree': 2})

# Create facet to cell connectivity required to determine boundary facets.
domain.topology.create_connectivity(domain.topology.dim, domain.topology.dim)
domain.topology.create_connectivity(domain.topology.dim, domain.topology.dim-1)
domain.topology.create_connectivity(domain.topology.dim-1, domain.topology.dim)

#  Define facet normal
n = ufl.FacetNormal(domain)

# %% [markdown]
# # Material parameters

# %%
# A function for constructing spatially varying (piecewise-constant) material parameters

# Volume domains: 
# Physical Volume("Free_vol", 44) = {1, 2};
# //+
# Physical Volume("Clamp_vol", 45) = {3};
# //+
# Physical Volume("Air_vol", 46) = {10};


# Need some extra infrastructure for the spatially-discontinuous material property fields
Vmat = functionspace(domain, ("DG", 0)) # create a DG0 function space on the domain

# A function to "safely" assign a value for properties, 
# accounting for syntax differences between FEniCSx v0.9.0 and v0.8.0. 
# This ensures that the code can be run using either version.
def safe_prop_assign(prop, i, val):
    try: # use the FEniCSx v0.9.0 syntax
        prop.x.petsc_vec.setValueLocal(i, val)
    except: # backup in case the user is running FEniCSx v0.8.0
        prop.vector.setValueLocal(i, val)

def mat(prop_val_cond, prop_val_air):

    # Define an empty "prop" material parameter function,
    # which lives on the DG0 function space.
    prop = Function(Vmat)
    
    # Now, actualy assign the desired values of the given material property to the new field.
    #
    coords = Vmat.tabulate_dof_coordinates()
    #
    # loop over the coordinates and assign the relevant material property, 
    # based on the local cell tag number.
    for i in range(coords.shape[0]):
        if cell_tags.values[i] == 46: 
            safe_prop_assign(prop, i, prop_val_air) 
        else:
            safe_prop_assign(prop, i, prop_val_cond) 
    
    return prop

# %%
# Volume domains: 
# Physical Volume("Free_vol", 44) = {1, 2};
# //+
# Physical Volume("Clamp_vol", 45) = {3};
# //+
# Physical Volume("Air_vol", 46) = {10};


# %%
# Elasticity parameters
Eyoung = 68.9e6 # Young's modulus, kPa --- aluminum
anu = 0.3 # Poisson's ratio
Gshear   = mat(Eyoung/(2.0*(1.0+anu)), 1.0e-3) # Shear modulus, kPa
Kbulk    = mat(Eyoung/(3.0*(1.0-2.0*anu)), 1.0e-5) # 101.0/1000) # Bulk modulus, kPa

# Mass density
rho =  mat(2.713e-6, 0.0) # 2.713e3 kg/m^3 = 2.713e-6 kg/mm^3 --- aluminum

# Magnetization parameters 
#
# Vacuum permeability
mu0 = Constant(domain, 1.256e-6*1e3) # Vacuum permeability,  mN / A^2
#
# material permeability (paramagnetic response)

chi = mat(0.0, 0.0)# unitless magnetic susceptibility (iron=5000, copper=0, air)
mu  = mu0*(1.0 + chi) # magnetic permeability

# Conductivity value
sigma = mat(25.321e-3, 1.0e-13)   # 1/(muOhm mm) (Al, air). Recall 1 MS/m = 1e-3 1/(muOhm mm)

#  flag for whether the domain is the air or not
TairInd   = mat(0, 1)

# Generalized-alpha method parameters
alpha   = Constant(domain, 0.0) # Set \alpha=0, since we observe no spurious pressure oscillations
gamma   = Constant(domain, PETSc.ScalarType(0.5+alpha))
beta    = Constant(domain, PETSc.ScalarType((gamma+0.5)**2/4.))


# %% [markdown]
# # Simulation time-control related params

# %%
# Simulation time control-related params
t    = 0.0         # start time (s)
#
dt = 0.000001 # initial time increment (s)
#
Ttot   = 0.15       # total simulation time (s) 

# Compiler variables for time step
dk = Constant(domain, PETSc.ScalarType(dt))
dk_old = Constant(domain, PETSc.ScalarType(dt))

# magnitude of b-fields
by_mag = 0.055 # Teslas
bx_mag = 0.5 # Teslas

# Decay time constant
tau = 6.6e-3 # 6.6 ms

def bx_ramp_hold(t):

    return bx_mag


def by_decay(t):

    by_  = by_mag*float(exp(-t/tau))

    return by_


# %% [markdown]
# # Function spaces

# %%

U2 = element("Lagrange", domain.basix_cell(), 2, shape=(3,))  # For displacement
U2a = element("Lagrange", domain.basix_cell(), 2, shape=(3,))  # For vector potential
P1 = element("Lagrange", domain.basix_cell(), 2)  # For phi
P0 = element("DG", domain.basix_cell(), 0)  # For p
#
TH = mixed_element([U2, U2a, P1, P0])     # Mixed element
ME = functionspace(domain, TH)    # Total mixed function space for all DOFs
#
V1 = functionspace(domain, P1) # Scalar function space.
V2 = functionspace(domain, U2) # (disp.) Vector function space
V2a = functionspace(domain, U2a) # (magnetic potential) Vector function space

# Define actual functions with the required DOFs
w = Function(ME)
u, a, phi, p = split(w)  # displacement u, magnetic vector potential a, electrostatic potential phi

# A copy of functions to store values in the previous step for time-stepping
w_old = Function(ME)
u_old, a_old, phi_old, p_old = split(w_old)   
#
# functions for "older" dofs from two steps ago
w_old2 = Function(ME)
u_old2, a_old2, phi_old2, p_old2 = split(w_old2)

# Define test functions    
u_test, a_test, phi_test, p_test = TestFunctions(ME)  

# Define trial functions needed for automatic differentiation
dw = TrialFunction(ME)  

# Functions for storing the velocity and acceleration at prev. step
v_old = Function(V2)
acc_old = Function(V2)


# %% [markdown]
# # Initial conditions
# 
# - The initial conditions for degrees of freedom $\mathbf{u}$, $\mathbf{a}$, and $\Phi$ are zero everywhere.
# - These are imposed automatically, since we have not specified any non-zero initial conditions.

a1_0_expr = by_mag*x[2]/2.0   + 0.0
a2_0_expr = 0.0               - bx_mag*x[2]/2.0
a3_0_expr = -by_mag*x[0]/2.0  + bx_mag*x[1]/2.0
a0_vec  = ufl.as_vector([a1_0_expr, a2_0_expr, a3_0_expr])
#
expr_a0    = Expression(a0_vec,V2a.element.interpolation_points())
w_old.sub(1).interpolate(expr_a0)
w_old2.sub(1).interpolate(expr_a0)

# %% [markdown]
# # Subroutines for kinematics and constitutive equations

# %%
#------------------------------------------------------------- 
# Utility subroutines
#-------------------------------------------------------------
 
# Subroutine for a "safer" sqrt() function which avoids a divide by zero 
# when automatically differentiated. 
def safe_sqrt(x):
    return sqrt(x + 1.0e-16)

# Plane strain deformation gradient 
def F_calc(u_):
    
    dim = len(u_)                # dimension of problem (3)
    
    Id = Identity(dim)          # 3D Identity tensor
    
    F = Id + grad(u_)            # 3D Deformation gradient

    return F


# Subroutine for calculating the (mechanical) Cauchy stress
def T_mech_calc(u_, p_):
    
    F   = F_calc(u_)
    
    J = det(F)
    
    Fbar = J**(-1/3)*F
    
    Bbar = Fbar*Fbar.T
    
    T_mech = (1/J)* Gshear * dev(Bbar) - p_ * Kbulk * Identity(3)
    
    return T_mech


def Piola_calc(u_, a_, p_):
    
    F = F_calc(u_)
    
    J = det(F)
    
    T_mech = T_mech_calc(u_, p_)

    Piola = J * T_mech * inv(F.T)
    
    return  Piola


#---------------------------------------------------------------------
# Subroutine for updating  acceleration using the Newmark beta method:
# a = 1/(2*beta)*((u - u0 - v0*dt)/(0.5*dt*dt) - (1-2*beta)*a0)
#---------------------------------------------------------------------
def update_a(u, u_old, v_old, a_old):
    return (u-u_old-dk*v_old)/beta/dk**2 - (1-2*beta)/2/beta*a_old

#---------------------------------------------------------------------
# Subroutine for updating  velocity using the Newmark beta method
# v = dt * ((1-gamma)*a0 + gamma*a) + v0
#---------------------------------------------------------------------
def update_v(a, u_old, v_old, a_old):
    return v_old + dk*((1-gamma)*a_old + gamma*a)

#---------------------------------------------------------------------
# alpha-method averaging function
#---------------------------------------------------------------------
def avg(x_old, x_new, alpha):
    return alpha*x_old + (1-alpha)*x_new


# %% [markdown]
# # Evaluate kinematics and constitutive relations

# %%
# Get acceleration and velocity at end of step
acc_new = update_a(u, u_old, v_old, acc_old)
v_new = update_v(acc_new, u_old, v_old, acc_old)

# get avg (u,p) fields for generalized-alpha method
u_avg   = avg(u_old, u, alpha)
a_avg   = avg(a_old, a, alpha)
phi_avg = avg(phi_old, phi, alpha)
p_avg   = avg(p_old, p, alpha)

# kinematical quantities
F  = F_calc(u_avg)
J  = det(F)
C = F.T*F

# EM potential-derived fields
#
# Magnetic flux density
bR = ufl.curl(a_avg)
#
# \dot{a}_R
adot = ((2.0*dk + dk_old)*dk_old*a_avg - ((dk+dk_old)**2)*a_old + (dk**2)*a_old2 ) / ((dk*dk_old)*(dk + dk_old))
#
# Electromotive intensity
phiNorm = Constant(domain, PETSc.ScalarType(1.0e6))
calE_R = - adot - phiNorm * ufl.grad(phi_avg)

# Eddy currents
je   = J*sigma*ufl.inv(C)*calE_R

# Magnetic field
hR =  1.0/(mu*J)*C*ufl.curl(a_avg)

# Mechanical Piola stress
Piola = Piola_calc(u_avg, a_avg, p_avg)

# Lorentz force
lorentz = inv(F.T)*ufl.cross(je, bR)

# %%
# Volume domains: 
# Physical Volume("Free_vol", 44) = {1, 2};
# //+
# Physical Volume("Clamp_vol", 45) = {3};
# //+
# Physical Volume("Air_vol", 46) = {10};

# domain of the conducting body \B
dx_eddy = dx(44) + dx(45)

# The weak form for the equilibrium equation
# 
Res_0  =  inner( Piola, grad(u_test))*dx  \
          + inner(rho * acc_new, u_test)*dx \
         - inner(lorentz, u_test)*dx_eddy

# Weak form of Ampere's Law
Res_1 = inner( hR, ufl.curl(a_test) )*dx \
       - inner( je, a_test )*dx_eddy \
       + inner(  1.0/(float(mu0)*J) * ufl.div(a_avg), ufl.div(a_test) ) * dx 

# weak form of charge balance law, in conductor
Res_2 = inner(   je, ufl.grad(phi_test)) * dx 

Res_3 = inner( (J-1) + p, p_test) * dx

Res = (1/float(Eyoung))*Res_0 + Res_1 + Res_2 + Res_3

# %% [markdown]
# # Set-up output files

# %%
# results file name
results_name = "TEAM12_3D_eddy_bending_bx0p5T"

# Function space for projection of results
P1 = element("DG", domain.basix_cell(), 1)
VV1 = fem.functionspace(domain, P1) # linear scalar function space
#
U1 = element("DG", domain.basix_cell(), 1, shape=(3,)) 
VV2 = fem.functionspace(domain, U1) # linear Vector function space
#
T1 = element("DG", domain.basix_cell(), 1, shape=(3,3)) 
VV3 = fem.functionspace(domain, T1) # linear tensor function space
#
U1a = element("DG", domain.basix_cell(), 1, shape=(3,)) 
VV2a = fem.functionspace(domain, U1a) # linear Vector function space

# For visualization purposes, we need to re-project the stress tensor onto a linear function space before 
# we write it (and its components and the von Mises stress, etc) to the VTX file. 
#
# This is because the stress is a complicated "mixed" function of the (quadratic Lagrangian) displacements
# and the (quadrature representation) plastic strain tensor and scalar equivalent plastic strain. 
#
# Create a linear problem for projecting the stress tensor onto the linear tensor function space VV3.

dx_proj = ufl.Measure('dx', domain=domain, subdomain_data=cell_tags, metadata={'quadrature_degree': 10})

def setup_projection(u, V):

    trial = ufl.TrialFunction(V)
    test  = ufl.TestFunction(V)  

   #  u_mod = conditional(gt(x[0], hmin), u, 0.0*u) 

    a_proj = ufl.inner(trial, test)*dx
    L_proj = ufl.inner(u, test)*dx_eddy

    # air_dofs = fem.locate_dofs_topological(V, cell_tags.dim, cell_tags.find(11))
    # bcs = [dirichletbc(0.0, air_dofs, V) ] 

    projection_problem = dolfinx.fem.petsc.LinearProblem(a_proj, L_proj, [], \
       petsc_options={"ksp_type": "cg", "ksp_rtol": 1e-16, "ksp_atol": 1e-16, "ksp_max_it": 1000})
    
    return projection_problem


# %%
# primary fields to write to output file
u_vis      = Function(VV2a, name="disp")
a_vis      = Function(VV2a, name="a_R")
phi_vis    = Function(VV1, name="Phi")
p_vis      = Function(VV1, name="p")

# %%
# Visualizing eddy currents
je_vis = Function(VV2a, name="je")
je_expr = Expression(F*je/J, VV2a.element.interpolation_points())

# Visualizing the power dissipation density
je_diss = dot(je, calE_R)
je_diss_vis = Function(VV1, name="je diss")
je_diss_expr = Expression(je_diss, VV1.element.interpolation_points())

# visualizing the magnetic b-field
b_vis = Function(VV2a, name="b_sp")
b_sp = F*ufl.curl(a)/J
b_expr = Expression(b_sp, VV2a.element.interpolation_points())

# Visualizing the (spatial) Lorentz force ( j x b ) term
je_sp = F*je/J
lorentz_sp = ufl.cross(je_sp, b_sp)
#
lorentz_proj = setup_projection(lorentz_sp, VV2)
lorentz_temp = lorentz_proj.solve()
#
lorentz_vis = Function(VV2, name="Lorentz")
lorentz_expr = Expression(lorentz_temp, VV2.element.interpolation_points())

# visualizing adot
adot_vis = Function(VV2a, name="adot")
adot_expr = Expression(adot, VV2a.element.interpolation_points())

# visualizing the air domain flag
air_flag_vis = Function(VV1, name="Air_flag")
air_flag_expr = Expression(TairInd, VV1.element.interpolation_points())
air_flag_vis.interpolate(air_flag_expr)

# Visualizing the shear modulus
Gshear_vis = Function(VV1, name="Gshear")
Gshear_expr = Expression(Gshear, VV1.element.interpolation_points())
Gshear_vis.interpolate(Gshear_expr)


# %%
# set up the output VTX files.
file_results = VTXWriter(
    MPI.COMM_WORLD,
    "results/" + results_name + ".bp",
    [  # put the functions here you wish to write to output
        u_vis,  a_vis, phi_vis, p_vis,  Gshear_vis, je_vis, je_diss_vis, b_vis, lorentz_vis, adot_vis, air_flag_vis, 
    ],
    engine="BP4",
)

def writeResults(t):
    
    # Update the output fields before writing to VTX.
    #
    u_vis.interpolate(w.sub(0))
    a_vis.interpolate(w.sub(1))
    phi_vis.interpolate(w.sub(2))
    p_vis.interpolate(w.sub(3))
    je_diss_vis.interpolate(je_diss_expr)
    je_vis.interpolate(je_expr)
    b_vis.interpolate(b_expr)
    lorentz_temp = lorentz_proj.solve()
    lorentz_vis.interpolate(lorentz_expr)
    adot_vis.interpolate(adot_expr)
       
    # Finally, write output fields to VTX.
    #
    file_results.write(t) 

# %% [markdown]
# # Infrastructure for pulling out time history data (force, displacement, etc.)

# %%
# infrastructure for evaluating functions at a certain point efficiently
pointForEval = np.array([487.0, 3.175, 100.0])
bb_tree = dolfinx.geometry.bb_tree(domain,domain.topology.dim)
cell_candidates = dolfinx.geometry.compute_collisions_points(bb_tree, pointForEval)
colliding_cells = dolfinx.geometry.compute_colliding_cells(domain, cell_candidates, pointForEval).array

# computing the total dissipation in the conductive domains
dissTot  = je_diss*dt*dx_eddy
dissForm = fem.form(dissTot) 

# Computing the total eddy current in the conductive part
# the current integration surface is 49.
flag = (1.0 + dot(n, ufl.as_vector([1,0,0])))/2.0
currentTot  = flag('+')*dot(je('+'), n('+'))*ds(49) + flag('-')*dot(je('-'), n('-'))*ds(49)
currentForm = fem.form(currentTot) 

# %% [markdown]
# ## Boundary condtions

# %%
# surface domains: 
# //+
# Physical Surface("Clamp_surfs", 47) = {17, 15, 14, 16, 1};
# //+
# Physical Surface("Air_surfs", 48) = {21, 18, 22, 23, 19, 20};
# //+
# Physical Surface("current_surf", 49) = {7};
# //+
# Physical Surface("maxwell_surf", 50) = {13, 2, 11, 10, 9, 3, 6, 4, 5};

# %%

# Find the specific DOFs which will be constrained.

airSurf_u_dofs = fem.locate_dofs_topological((ME.sub(0), V2), facet_tags.dim, facet_tags.find(48))
airSurf_a_dofs = fem.locate_dofs_topological((ME.sub(1), V2a), facet_tags.dim, facet_tags.find(48))
airSurf_phi_dofs = fem.locate_dofs_topological((ME.sub(2), V1), facet_tags.dim, facet_tags.find(48))
#
clampSurf_u_dofs = fem.locate_dofs_topological((ME.sub(0), V2), facet_tags.dim, facet_tags.find(47))

# Zero bc quantity for u
Vu_0, submap = ME.sub(0).collapse()
zero_u   = Function(Vu_0)
zero_u.interpolate(lambda x: np.stack(( np.zeros(x.shape[1]), np.zeros(x.shape[1]),  np.zeros(x.shape[1]) ) ) )

# Zero bc quantity for phi
Vphi_0, submap = ME.sub(2).collapse()
zero_phi   = Function(Vphi_0)
zero_phi.interpolate(lambda x:  np.zeros(x.shape[1])) 

# Time-varying constant for the ramped up magnetic vector potential component
by_applied = Constant(domain, PETSc.ScalarType(by_decay(t)))
bx_applied = Constant(domain, PETSc.ScalarType(bx_ramp_hold(t)))

# Expressions for magnetic vector potential BCs
a1_expr = by_applied*x[2]/2.0   + 0.0
a2_expr = 0.0                   - bx_applied*x[2]/2.0
a3_expr = -by_applied*x[0]/2.0  + bx_applied*x[1]/2.0
a_vec  = ufl.as_vector([a1_expr, a2_expr, a3_expr])
#
expr_a  = Expression(a_vec,V2.element.interpolation_points())
func_a  = Function(V2)
func_a.interpolate(expr_a)

# building Dirichlet BCs
bcs_0 = dirichletbc(zero_u, airSurf_u_dofs, ME.sub(0)) # displacement fixed, air domain surface
bcs_1 = dirichletbc(func_a, airSurf_a_dofs, ME.sub(1)) # magnetic vector potential expression, air domain surface
bcs_2 = dirichletbc(zero_phi, airSurf_phi_dofs, ME.sub(2)) # phi fixed, air domain surface
#
bcs_3 = dirichletbc(zero_u, clampSurf_u_dofs, ME.sub(0)) # displacement fixed, clamped surface


bcs = [bcs_0, bcs_1, bcs_2, bcs_3]


# %% [markdown]
# ## Define the nonlinear variational problem

# %%
# a function for re-building the solver when the j_s distribution is changed.
def build_problem():

    # Automatic differentiation tangent:
    J_tan = derivative(Res, w, dw)

    # Set up nonlinear problem
    problem = NonlinearProblem(Res, w, bcs, J_tan)
    
    # the global newton solver and params
    solver = NewtonSolver(MPI.COMM_WORLD, problem)
    solver.convergence_criterion = "incremental"
    solver.rtol = 1e-6
    solver.atol = 1e-6
    solver.max_it = 25
    solver.error_on_nonconvergence = False

    return solver


# %% [markdown]
# ##  Start calculation loop

# %%
# Give the step a descriptive name
step = "Actuate"

# Variables for storing time history
totSteps = 1000000
timeHist0 = np.zeros(shape=[totSteps])
timeHist1 = np.zeros(shape=[totSteps]) 
timeHist2 = np.zeros(shape=[totSteps])  
timeHist3 = np.zeros(shape=[totSteps]) 
timeHist4 = np.zeros(shape=[totSteps]) 

#Iinitialize a counter for reporting data
ii=0

#  Set up temporary "helper" functions and expressions 
#  for updating the velocity and acceleration.
v_temp = Function(V2)
acc_temp = Function(V2)
#
v_expr = Expression(v_new,V2.element.interpolation_points())
acc_expr = Expression(acc_new,V2.element.interpolation_points())


# Write initial state to file
writeResults(t=0.0)   

# Build the problem and associated solver
solver = build_problem()

# print a message for simulation startup
mprint("------------------------------------")
mprint("Simulation Start")
mprint("------------------------------------")
# Store start time 
startTime = datetime.now()

# Time-stepping solution procedure loop
while (round(t, 9) <= Ttot):
     
    # increment time
    t += float(dk) 
    
    # update time variables in time-dependent BCs 
    #
    # Update the applied b-field value
    t_avg = t - float(alpha*dk)
    by_applied.value = by_decay(t_avg)
    bx_applied.value = bx_ramp_hold(t_avg)
    
    # re-interpolate the magnetic vector potential
    func_a.interpolate(expr_a)
    
    # Solve the problem
    (iter, converged) = solver.solve(w) 
    
    # Collect results from MPI ghost processes
    w.x.scatter_forward()
        
    # Adaptive time-stepping procedure: 
    #
    # Check if solver converges
    if converged: 
        
        # increment counter
        ii += 1   
        
        # reset non-convergence counter
        nonCon = 0
        
        # Print progress of calculation periodically
        if ii%1 == 0:      
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            mprint("Step: {} |   Increment: {} | Iterations: {}".format(step, ii, iter))
            mprint("dt: {} s | Simulation Time: {} s  of  {} s".format(round(dt, 4), round(t,4), Ttot))
            mprint()   
            
        # Write output to file
        writeResults(t)
        
        # Store time history variables at this time  

        # Store time history variables
        if len(colliding_cells)>=1:
            timeHist0[ii] = t # current time 
            #
            timeHist1[ii] = float(by_applied.value) # time history of applied b-field
            #
            timeHist4[ii] = w.sub(0).sub(1).eval([487.0, 3.175, 100.0],colliding_cells[0])[0] # time history of b-field
            #
            plotFlag = True
        else:
            plotFlag = False
        #
        timeHist2[ii] = domain.comm.allreduce(fem.assemble_scalar(currentForm), op=MPI.SUM)
        # time integral of dissipation in the channel
        timeHist3[ii] = timeHist3[ii-1] + domain.comm.allreduce(fem.assemble_scalar(dissForm), op=MPI.SUM)  

        # Save the displacement and current data to a .csv file each time step
        # (but only on one process)
        if plotFlag:
            histData = np.array([timeHist0[0:ii], timeHist4[0:ii], timeHist2[0:ii]])
            np.savetxt("results/bx0p5_data.csv", histData, delimiter=",")
        
        # update DOFs and internal variables 
        #
        # interpolate the values of the internal variables into their "temp" functions
        v_temp.interpolate(v_expr)
        acc_temp.interpolate(acc_expr)
        #
        # Update DOFs for next step
        w_old2.x.array[:] = w_old.x.array
        w_old.x.array[:] = w.x.array
        #
        v_old.x.array[:] = v_temp.x.array[:]
        acc_old.x.array[:] = acc_temp.x.array[:]
            
        # store the time step size
        dk_old.value = dk.value
        
        # Iteration-based adaptive time-stepping
        dtMax = Ttot/100
        if (iter<=3):
            dt = 1.5*dt
            dk.value = dt
        elif iter>=5:
            dt = 0.5*dt
            dk.value = dt
        
        # limit maximum step size
        if dt>dtMax:
            dt = dtMax
            dk.value = dt
            
    # If solver does not  converge, then do not  save results and try a smaller dt 
    else:
        
        # increment non-convergence counter
        nonCon += 1
        
        # if >5 consecutive non-convergences, kill the simulation.
        if nonCon>5:
            mprint("Too many non-converging increments. Killing simulation.")
            break
                # back up in time
        t = t - float(dk)
        
        # cut back on dt
        dt = dt/10
        dk.value = dt
        
        # Reset DOFs for next step
        w.x.array[:] = w_old.x.array 
        
        # re-build the problem
        solver = build_problem()
        
        # print out the new dt value
        mprint("Attempting dt = {}".format(dt))
        mprint()
            
            
# close the output file.
file_results.close()
         
# End analysis
mprint("-----------------------------------------")
mprint("End computation")                 
# Report elapsed real time for the analysis
endTime = datetime.now()
elapseTime = endTime - startTime
mprint("------------------------------------------")
mprint("Elapsed real time:  {}".format(elapseTime))
mprint("------------------------------------------")


# %% [markdown]
# # Plot results

# %%
if plotFlag: 

    # Set plot font to size 14
    font = {'size'   : 14}
    plt.rc('font', **font)

    # Get array of default plot colors
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors     = prop_cycle.by_key()['color']

    # Only plot as far as we have time history data
    ind = np.argmax(timeHist0) +1

    plt.figure()
    plt.plot(timeHist0[0:ind], timeHist2[0:ind], linewidth=2.0, color=colors[1]) #, marker='.')
    plt.axis('tight')
    plt.xlabel(r"Time (s)")
    plt.ylabel(r"Eddy current (A)")
    plt.grid(linestyle="--", linewidth=0.5, color='b')

    fig = plt.gcf()
    fig.set_size_inches(7,5)
    plt.tight_layout()
    plt.savefig("results/eddy_current_curve_MPI_bx0p5.png", dpi=600)



    plt.figure()
    plt.plot(timeHist0[0:ind], -timeHist4[0:ind], linewidth=2.0, color=colors[3]) #, marker='.')
    plt.axis('tight')
    plt.xlabel(r"Time (s)")
    plt.ylabel(r"Displacement (mm)")
    plt.grid(linestyle="--", linewidth=0.5, color='b')

    fig = plt.gcf()
    fig.set_size_inches(7,5)
    plt.tight_layout()
    plt.savefig("results/tip_disp_curve_MPI_bx0p5.png", dpi=600)

