using WAVI 
function driver()

#
#Grid and boundary conditions
#
nx = 150
ny = 25
nσ = 4
x0 = 0.0
y0 = -25000.0
dx = 2000.0
dy = 2000.0
h_mask=trues(nx,ny)

#
# boundary conditions
#
u_iszero = falses(nx+1,ny); 
u_iszero[1,:].=true #xero u flow at x = 0
v_iszero=falses(nx,ny+1); 
v_iszero[:,1].=true;        #zero v flow at y =  25km 
v_iszero[:,end].=true       #zero v flow at y = -25km
u_iszero[:,1].=true;
u_iszero[:,end].=true;      #no u velocity at lateral boundaries

grid = Grid(nx = nx, 
            ny = ny,   
            nσ = nσ, 
            x0 = x0, 
            y0 = y0, 
            dx = dx, 
            dy = dy,
            h_mask = h_mask, 
            u_iszero = u_iszero, 
            v_iszero = v_iszero)

#
#Bed 
#
bathy=Array{Float64}(undef,nx,ny);
read!("bathy.bin",bathy)
bathy.=ntoh.(bathy)

# 
# initial conditions
#
h=Array{Float64}(undef,nx,ny);
read!("h_init.bin",h)
h.=ntoh.(h)
initial_conditions = InitialConditions(initial_thickness = h)

#
#solver parameters
#
maxiter_picard = 1
solver_params = SolverParams(maxiter_picard = maxiter_picard)

#
#Physical parameters
#
accumulation =Array{Float64}(undef,nx,ny);
read!("accumulation.bin",accumulation)
accumulation.=ntoh.(accumulation)

weertman_c = 368.0

glen_a_ref =Array{Float64}(undef,nx,ny);
read!("glen_a_ref.bin",glen_a_ref)
glen_a_ref.=ntoh.(glen_a_ref)

params = Params(accumulation_rate = accumulation,
				  glen_a_ref = glen_a_ref,
				  weertman_c = weertman_c)


#
# Melt rate
#
# time_fpath =  joinpath(dirname(@__FILE__), "ref_time.bin")
# depth_fpath =  joinpath(dirname(@__FILE__), "ref_depth.bin")
# Ta_fpath =  joinpath(dirname(@__FILE__), "forcing_ambient_temperature.bin")
# Sa_fpath =  joinpath(dirname(@__FILE__), "forcing_ambient_salinity.bin")

#load files and put in correct format
# nt = convert(Int, floor(filesize(time_fpath)/8)) #8 bit binary
# ref_time = Array{Float64}(undef, nt)
# read!(time_fpath, ref_time)
# ref_time = ntoh.(ref_time)

# nz = convert(Int, floor(filesize(depth_fpath)/8)) #8 bit binary
# ref_depth = Array{Float64}(undef, nz)
# read!(depth_fpath, ref_depth)
# ref_depth = ntoh.(ref_depth)

# Ta = Array{Float64}(undef, nt, nz)
# read!(Ta_fpath, Ta)
# Ta = ntoh.(Ta);

# Sa = Array{Float64}(undef, nt, nz)
# read!(Sa_fpath, Sa)
# Sa = ntoh.(Sa);
# melt = QuadraticForcedMeltRate(ref_depth = ref_depth,
#                            ref_time = ref_time,
#                            Ta = Ta,
#                            Sa = Sa,
#                            γT = 1.1*1e-3,
#			    melt_partial_cell = true);

#
#make the model
#
model = Model(grid = grid,
              bed_elevation = bathy, 
              params = params, 
              solver_params = solver_params,
	      initial_conditions = initial_conditions)#,
	     # melt_rate = melt);

#
#timestepping parameters
#
niter0 = 0
dt = 0.01
end_time = 500.
chkpt_freq = 10.
pchkpt_freq = 10.
timestepping_params = TimesteppingParams(niter0 = niter0, 
                                         dt = dt, 
                                         end_time = end_time, 
                                         chkpt_freq = chkpt_freq, 
                                          pchkpt_freq = pchkpt_freq)

#
#output parameters
#
outputs = (h = model.fields.gh.h,
           u = model.fields.gh.u,
           v = model.fields.gh.v,
           b = model.fields.gh.b,
           s = model.fields.gh.s,
           grfrac = model.fields.gh.grounded_fraction,
           m = model.fields.gh.basal_melt)

output_freq = 10.
output_params = OutputParams(outputs = outputs, 
                            output_freq = output_freq,
                            output_format = "mat",
                            dump_vel = true,
                            zip_format = "nc",
                            output_start = true)

#
# assemble the simulation
#
simulation = Simulation(model = model, 
                        timestepping_params = timestepping_params, 
                        output_params = output_params)
                
#
#perform the simulation
#
run_simulation!(simulation)

return simulation
end

driver()

