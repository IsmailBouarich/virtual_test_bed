# ==============================================================================
# High Temperature Transient Facility (HTTF)
# Main App solve of the HTTF steady state 
# ------------------------------------------------------------------------------
# Idaho Falls, INL, 08/2024
# Author(s): Ismail Bouarich, Dr. Lise Charlot
# ==============================================================================
# This input file handles 3-D porous medium model and communicates back and forth with subscale models

!include common.i

# physics parameters for the porous medium model 
inlet_temperature = 500 #K
outlet_pressure = 7e5 # Pa
power =2.2e6 # W
pin_power_density = '${fparse power/210/(A_hexa * heated_length)}' #power density in a single fuel rod pin

#Solid material properties 
rho_s = 1870.006372
cp_s = 1711.8617
k_s = 82.5780914

#Values of porosity for each region
por_0 = 0.00000000
por_1 = 0.04047856
por_2 = 0.02023928
por_3 = 0.03822975
por_4 = 0.16866067
por_5 = 0.14842139
por_6 = 0.07421070
por_7 = 0.02023928
por_8 = 0.00000000
por_9 = 0.02811011

#Air to volume fraction for each region
a_to_v_0 = 0.00000000
a_to_v_1 = 8.49943548
a_to_v_2 = 8.49943548
a_to_v_3 = 14.16572581
a_to_v_4 = 42.49717742
a_to_v_5 = 39.66403226
a_to_v_6 = 25.49830645
a_to_v_7 = 8.49943548
a_to_v_8 = 0.00000000
a_to_v_9 = 7.08286290

#Hydraulic diameter for each region
Dh_0 = 0.00000000
Dh_1 = 0.01905000
Dh_2 = 0.00952500
Dh_3 = 0.01079500
Dh_4 = 0.01587500
Dh_5 = 0.01496786
Dh_6 = 0.01164167
Dh_7 = 0.00952500
Dh_8 = 0.00000000
Dh_9 = 0.01587500

#Block categorization
flow_blocks = 'core_iref_cring core_iref_oring_cen core_iref_oring_cor
               core_fuel_inner core_fuel_outer
               core_oref_iring_cor core_oref_iring_con core_oref_oring

               br_iref_cring br_iref_oring_cen br_iref_oring_cor
               br_fuel_inner br_fuel_outer
               br_oref_iring_cor br_oref_iring_con br_oref_oring

               tr_iref_cring tr_iref_oring_cen tr_iref_oring_cor
               tr_fuel_inner tr_fuel_outer
               tr_oref_iring_cor tr_oref_iring_con tr_oref_oring

               helium_gap inlet_plenum'
flow_blocks_core = 'core_iref_cring core_iref_oring_cen core_iref_oring_cor
               core_fuel_inner core_fuel_outer
               core_oref_iring_cor core_oref_iring_con core_oref_oring

               br_iref_cring br_iref_oring_cen br_iref_oring_cor
               br_fuel_inner br_fuel_outer
               br_oref_iring_cor br_oref_iring_con br_oref_oring

               tr_iref_cring tr_iref_oring_cen tr_iref_oring_cor
               tr_fuel_inner tr_fuel_outer
               tr_oref_iring_cor tr_oref_iring_con tr_oref_oring

                inlet_plenum'

solid_blocks = 'core_iref_iring core_oref_cring
                br_iref_iring   br_oref_cring
                tr_iref_iring   tr_oref_cring
                permanent_reflector core_barrel_cell'

all_blocks = 'core_iref_cring core_iref_oring_cen core_iref_oring_cor
              core_fuel_inner core_fuel_outer
              core_oref_iring_cor core_oref_iring_con core_oref_oring

              br_iref_cring br_iref_oring_cen br_iref_oring_cor
              br_fuel_inner br_fuel_outer
              br_oref_iring_cor br_oref_iring_con br_oref_oring

              tr_iref_cring tr_iref_oring_cen tr_iref_oring_cor
              tr_fuel_inner tr_fuel_outer
              tr_oref_iring_cor tr_oref_iring_con tr_oref_oring

              helium_gap inlet_plenum

              core_iref_iring core_oref_cring
              br_iref_iring   br_oref_cring
              tr_iref_iring   tr_oref_cring
              permanent_reflector core_barrel_cell'

power_blocks = 'core_iref_oring_cen core_iref_oring_cor
                core_fuel_inner core_fuel_outer
                core_oref_iring_cor core_oref_iring_con'

# Creates a mesh from the given file                
[Mesh]
  [fmesh]
    type = FileMeshGenerator
    file = 'HTTF_mesh_complete_in.e'
  []
  # file = DCC_phase_1.e
  # distribution = serial 
[]

[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
  allow_initial_conditions_with_restart = true
[]

[GlobalParams]
  # rhie_chow_user_object = pins_rhie_chow_interpolator
[]

[Debug]
  # show_material_props = true
  # show_var_residual_norms = true
[]
# Solves the Navier-Stokes equations 
[Modules]
  [NavierStokesFV]
    # general input parameters
    compressibility = weakly-compressible
    block = ${flow_blocks}
    gravity = '0 0 -9.81'

    # porous media parameters
    porous_medium_treatment = true
    porosity = porosity
    friction_types = 'darcy forchheimer'
    friction_coeffs = 'Darcy_coefficient Forchheimer_coefficient'
    standard_friction_formulation = true


    # discretization parameters
    pressure_face_interpolation = average
    momentum_advection_interpolation = upwind
    mass_advection_interpolation = upwind

    # material properties
    density = 'rho'
    dynamic_viscosity = 'mu'

    # initial conditions
    initial_velocity = '0.0 0.0 0'
    initial_pressure = '${outlet_pressure}'
    # velocity_variable = 'superficial_vel_x superficial_vel_y superficial_vel_z'
    # pressure_variable = 'pressure'

    # boundary conditions
    inlet_boundaries = 'inlet_top'
    # momentum_inlet_types = 'flux-mass'
    # flux_inlet_pps = 'mdot_in'
    # flux_inlet_directions = '0 0 -1'
    momentum_inlet_types = 'fixed-velocity'
    momentum_inlet_function = '0.0 0.0 vel_in_fn'

    #Wall parameters
    wall_boundaries = 'wall'
    momentum_wall_types = 'slip'

    outlet_boundaries = ' iref_cring_bot
                          iref_oring_bot
                          iref_oring_cor_bot
                          fuel_inner_bot
                          fuel_outer_bot
                          oref_iring_cor_bot
                          oref_iring_con_bot
                          oref_oring_bot'
    momentum_outlet_types = 'fixed-pressure
                             fixed-pressure
                             fixed-pressure
                             fixed-pressure
                             fixed-pressure
                             fixed-pressure
                             fixed-pressure
                             fixed-pressure'
    pressure_function = 'pressure_fn
                         pressure_fn
                         pressure_fn
                         pressure_fn
                         pressure_fn
                         pressure_fn
                         pressure_fn
                         pressure_fn'
  []
[]
#Variables solved in the main app
[Variables]
  #Fluid temperature
  [T_fluid]
    type = INSFVEnergyVariable
    block = '${flow_blocks}'
    initial_condition = ${inlet_temperature}
  []
  #Solid temperature 
  [T_solid]
    type = INSFVEnergyVariable
    block = '${all_blocks}'
    initial_condition = ${inlet_temperature}
  []
[]
#Finite volume kernels used to solve for 3-D heat conduction 
[FVKernels]
  #This kernel implements time derivative term of the fluid heat conduction equation
  [temp_fluid_time_core]
    type = PINSFVEnergyTimeDerivative
    variable = T_fluid
    rho = 'rho'
    cp = 'cp'
    is_solid = false
    porosity = 'porosity'
    block = '${flow_blocks}'
  []
  #This kernel implements an advection term for the fluid heat transport equation
  [temp_fluid_advection_core]
    type = PINSFVEnergyAdvection
    variable = T_fluid
    rhie_chow_user_object = pins_rhie_chow_interpolator
    advected_interp_method = skewness-corrected
    block = '${flow_blocks}'
  []
  #This kernel implements a diffusion term for the fluid energy equation
  [temp_fluid_conduction_core]
    type = PINSFVEnergyAnisotropicDiffusion
    variable = 'T_fluid'
    kappa = 'kappa'
    kappa_interp_method = 'harmonic'
    porosity = 'porosity'
    block = ${flow_blocks}
  []
  #This kernel implements a volumetric convection heat transfer term for the fluid and solid energy equation converting the effect of solid temperature on fluid temperature
  [temp_solid_to_fluid]
    type = PINSFVEnergyAmbientConvection
    variable = 'T_fluid'
    T_fluid = 'T_fluid'
    T_solid = 'T_solid'
    is_solid = false
    h_solid_fluid = 'htc_vol'
    block = '${flow_blocks}'
  []
  #This kernel implements time derivative term of the solid heat conduction equation
  [temp_solid_time_core]
    type = PINSFVEnergyTimeDerivative
    variable = T_solid
    rho = 'rho_s'
    cp = 'cp_s'
    is_solid = true
    porosity = 'porosity'
    block = '${all_blocks}'
  []
  #This kernel implements a diffusion term for the solid energy equation
  [temp_solid_conduction_core]
    type = PINSFVEnergyAnisotropicDiffusion
    variable = 'T_solid'
    kappa = 'kappa_s'
    kappa_interp_method = 'harmonic'
    porosity = 'porosity'
    block = '${all_blocks}'
  []
  #This kernel implements the power of the heater rods
  [temp_solid_source]
    type = FVCoupledForce
    variable = T_solid
    v = 'heat_source'
    block = '${all_blocks}'
  []
  #This kernel implements a volumetric convection heat transfer term for the fluid and solid energy equation converting the effect of fluid temperature on solid temperature
  [temp_fluid_to_solid]
    type = PINSFVEnergyAmbientConvection
    variable = 'T_solid'
    T_fluid = 'T_fluid'
    T_solid = 'T_solid'
    is_solid = true
    h_solid_fluid = 'htc_vol'
    block = '${flow_blocks}'
  []
[]
#Finite volume boundary conditions for 3-D heat conduction in the porous medium model
[FVBCs]
  #Inlet fluid temperature boundary condition
  [T_fluid_in]
    type = FVFunctionDirichletBC
    variable = 'T_fluid'
    function = '${inlet_temperature}'
    boundary = 'inlet_top'
  []
  #Outlet solid temperature boundary condition
  [T_solid_out]
    type = FVFunctorConvectiveHeatFluxBC
    variable = 'T_solid'
    heat_transfer_coefficient = 10.0
    T_bulk = 'T_fluid'
    T_solid = 'T_solid'
    is_solid = true
    boundary = 'iref_cring_bot
                iref_oring_bot
                iref_oring_cor_bot
                fuel_inner_bot
                fuel_outer_bot
                oref_iring_cor_bot
                oref_iring_con_bot
                oref_oring_bot'
  []
  #Core barrel solid temperature boundary condition
  [T_solid_outer_barrel]
    type = FVFunctorConvectiveHeatFluxBC
    variable = 'T_solid'
    heat_transfer_coefficient = 10.0
    T_bulk = ${inlet_temperature}
    T_solid = 'T_solid'
    is_solid = true
    boundary = 'outer_core_barrel'
  []
[]
#Necessary functions for the system
[Functions]
  #Inlet velocity function
  [vel_in_fn]
    type = PiecewiseLinear
    x = '0 30'
    y = '0 -0.823' #0.823 originally
  []
  #Pressure function
  [pressure_fn]
    type = ConstantFunction
    value = ${outlet_pressure}
  []
  #Mass flow rate function
  [mdot_fn]
    type = ConstantFunction
    value = 1
  []
  #Power function
  [heater_power]
    type = ConstantFunction
    value = ${power}
  []
  #Power density in a hexagonal region
  [power_density_rod]
    type = ParsedFunction
    expression = 'power/(210 * A_hexa * heated_length)'
    symbol_names = 'power A_hexa heated_length'
    symbol_values = 'heater_power ${A_hexa} 1.982'
  []
[]  
#Implements heat conduction properties of the solid materials in the system
[SolidProperties]
  #Graphite thermal properties
  [G348_thermal]
    type = ThermalFunctionSolidProperties
    k = k_G348
    cp = cp_G348
    rho = rho_G348
  []
  #Greencast thermal properties
  [greencast_thermal]
    type = ThermalFunctionSolidProperties
    k = k_greencast_fn
    cp = cp_greencast_fn
    rho = rho_greencast
  []
[]
#Implements helium properties in the system
[FluidProperties]
  [helium_obj]
    type = HeliumFluidProperties
  []
[]
#Implements material properties in the system
[FunctorMaterials]
  #Creates material properties (specific heat, thermal conductivity...l) from fluid properties
  [fluid_props_to_mat_props]
    type = GeneralFunctorFluidProps
    fp = helium_obj
    T_fluid = 'T_fluid'
    characteristic_length = 'characteristic_length'
    porosity = 'porosity'
    pressure = 'pressure'
    speed = 'interstitial_vel'
  []
  #Computes enthalpy
  [ins_fv]
    type = INSFVEnthalpyFunctorMaterial
    rho = 'rho'
    cp = 'cp'
    temperature = 'T_fluid'
  []
  #Implements material properties
  [sp_mat]
    type = ThermalSolidPropertiesFunctorMaterial
    temperature = T_solid
    sp = greencast_thermal
    specific_heat = cp_s
    density = rho_s
    thermal_conductivity = k_s
  []
  #Liquid thermal conductivity vector
  [vector_conductivity_liquid]
    type = ADGenericVectorFunctorMaterial
    prop_names = 'kappa'
    prop_values = '0 0 k'
  []
  #Transversal solid thermal conductivity
  [transversal_solid_conductivity]
    type = NSFVMixtureFunctorMaterial
    prop_names = 'k_t'
    phase_1_names = 'k'
    phase_2_names = 'k_s'
    phase_1_fraction = 'porosity'
    averaging_method = reciprocal-average
  []
  #Solid thermal conductivity vector
  [vector_conductivity_solid]
    type = ADGenericVectorFunctorMaterial
    prop_names = 'kappa_s'
    prop_values = 'k_t k_t k_s'
  []

  ## Porosity
  [porosity_functor_mat]
    type = ADPiecewiseByBlockFunctorMaterial
    prop_name = 'porosity'
    subdomain_to_prop_value = ' core_iref_cring 	${por_1}
                                core_iref_oring_cen 	${por_2}
                                core_iref_oring_cor	${por_3}
                                core_fuel_inner 	${por_4}
                                core_fuel_outer	${por_5}
                                core_oref_iring_cor 	${por_6}
                                core_oref_iring_con	${por_7}
                                core_oref_oring	${por_9}

                                br_iref_cring 	${por_1}
                                br_iref_oring_cen 	${por_2}
                                br_iref_oring_cor	${por_3}
                                br_fuel_inner 	${por_4}
                                br_fuel_outer ${por_5}
                                br_oref_iring_cor 	${por_6}
                                br_oref_iring_con ${por_7}
                                br_oref_oring	${por_9}

                                tr_iref_cring 	${por_1}
                                tr_iref_oring_cen 	${por_2}
                                tr_iref_oring_cor	${por_3}
                                tr_fuel_inner 	${por_4}
                                tr_fuel_outer ${por_5}
                                tr_oref_iring_cor 	${por_6}
                                tr_oref_iring_con ${por_7}
                                tr_oref_oring	${por_9}

                                helium_gap 	1.0
                                inlet_plenum	1.0

                                core_iref_iring  0.0
                                core_oref_cring  0.0
                                br_iref_iring  0.0
                                br_oref_cring  0.0
                                tr_iref_iring  0.0
                                tr_oref_cring  0.0
                                permanent_reflector  0.0
                                core_barrel_cell  0.0'
    block = ${all_blocks}
  []

  ## set characteristic length on each block
  [characteristic_length]
    type = PiecewiseByBlockFunctorMaterial
    prop_name = 'characteristic_length'
    subdomain_to_prop_value = ' core_iref_cring 	${Dh_1}
                                core_iref_oring_cen 	${Dh_2}
                                core_iref_oring_cor	${Dh_3}
                                core_fuel_inner 	${Dh_4}
                                core_fuel_outer	${Dh_5}
                                core_oref_iring_cor 	${Dh_6}
                                core_oref_iring_con	${Dh_7}
                                core_oref_oring	${Dh_9}

                                br_iref_cring 	${Dh_1}
                                br_iref_oring_cen 	${Dh_2}
                                br_iref_oring_cor	${Dh_3}
                                br_fuel_inner 	${Dh_4}
                                br_fuel_outer	${Dh_5}
                                br_oref_iring_cor 	${Dh_6}
                                br_oref_iring_con	${Dh_7}
                                br_oref_oring ${Dh_9}

                                tr_iref_cring 	${Dh_1}
                                tr_iref_oring_cen 	${Dh_2}
                                tr_iref_oring_cor	${Dh_3}
                                tr_fuel_inner 	${Dh_4}
                                tr_fuel_outer	${Dh_5}
                                tr_oref_iring_cor 	${Dh_6}
                                tr_oref_iring_con	${Dh_7}
                                tr_oref_oring	${Dh_9}

                                helium_gap 	0.01
                                inlet_plenum	1.504'
    block = ${flow_blocks}
  []
  #Material providing helium drag coefficients computed from Churchill correlations as functor
  [fuel_drag]
    type = FunctorChurchillDragCoefficients
    roughness = 0.0
    block = ${flow_blocks_core}
    multipliers = '10000 10000 1'
  []
  [f_he_gap]
     type = FunctorChurchillDragCoefficients
    roughness = 10.0
    block = helium_gap
  []
  #Populates each block with the correct power density
  [heat_source_fraction]
    type = PiecewiseByBlockFunctorMaterial
    prop_name = 'heat_source_fraction'
    subdomain_to_prop_value = ' core_iref_cring 	0
                                core_iref_oring_cen 	0.333333333333333
                                core_iref_oring_cor	0.666666666666667
                                core_fuel_inner 	3
                                core_fuel_outer	2.666666666666
                                core_oref_iring_cor 	1.33333333333333
                                core_oref_iring_con	0.333333333333333
                                core_oref_oring	0

                                br_iref_cring 	0
                                br_iref_oring_cen 	0
                                br_iref_oring_cor	0
                                br_fuel_inner 	0
                                br_fuel_outer	0
                                br_oref_iring_cor 	0
                                br_oref_iring_con	0
                                br_oref_oring	0

                                tr_iref_cring 	0
                                tr_iref_oring_cen 	0
                                tr_iref_oring_cor	0
                                tr_fuel_inner 	0
                                tr_fuel_outer	0
                                tr_oref_iring_cor 	0
                                tr_oref_iring_con	0
                                tr_oref_oring	0

                                helium_gap 	0
                                inlet_plenum	0

                                core_iref_iring	0
                                core_oref_cring	0
                                br_iref_iring	0
                                br_oref_cring	0
                                tr_iref_iring	0
                                tr_oref_cring	0
                                permanent_reflector	0
                                core_barrel_cell	0'
    block = ${all_blocks}
  []
  [heat_source]
    type = ParsedFunctorMaterial
    property_name = heat_source
    expression = 'power_density_rod * heat_source_fraction'
    functor_names = 'power_density_rod heat_source_fraction'
  []
  #Populates each block with the correct air to volume number
  [a_to_v]
    type = PiecewiseByBlockFunctorMaterial
    prop_name = 'a_to_v'
    subdomain_to_prop_value = ' core_iref_cring 	${a_to_v_1}
                                core_iref_oring_cen 	${por_2}
                                core_iref_oring_cor	${a_to_v_3}
                                core_fuel_inner 	${a_to_v_4}
                                core_fuel_outer	${a_to_v_5}
                                core_oref_iring_cor 	${a_to_v_6}
                                core_oref_iring_con	${a_to_v_7}
                                core_oref_oring	${a_to_v_9}

                                br_iref_cring 	${a_to_v_1}
                                br_iref_oring_cen 	${a_to_v_2}
                                br_iref_oring_cor	${a_to_v_3}
                                br_fuel_inner 	${a_to_v_4}
                                br_fuel_outer ${a_to_v_5}
                                br_oref_iring_cor 	${a_to_v_6}
                                br_oref_iring_con ${a_to_v_7}
                                br_oref_oring	${a_to_v_9}

                                tr_iref_cring 	${a_to_v_1}
                                tr_iref_oring_cen 	${a_to_v_2}
                                tr_iref_oring_cor	${a_to_v_3}
                                tr_fuel_inner 	${a_to_v_4}
                                tr_fuel_outer ${a_to_v_5}
                                tr_oref_iring_cor 	${a_to_v_6}
                                tr_oref_iring_con ${a_to_v_7}
                                tr_oref_oring	${a_to_v_9}

                                helium_gap 	12631.158439016303
                                inlet_plenum	1'
    block = ${all_blocks}
  []
  #Computes wall heat transfer coefficients using the Dittus-Boelter correlation
  [htc]
    type = FunctorDittusBoelterWallHTC

    block = ${flow_blocks}
  []
  #Converts the wall transfer coefficients into volume wall transfer coefficients
  [htc_vol]
    type = ADParsedFunctorMaterial
    functor_names = 'a_to_v wall_htc'
    property_name = 'htc_vol'
    expression = 'a_to_v * wall_htc'
  []

[]
#Creates relevant AuxVariables to the problem
[AuxVariables]
  #Power density variable
  [power_source]
    type = MooseVariableFVReal
    block = '${power_blocks}'
  []
  #Interstitial velocity variable
  [interstitial_vel]
    type = MooseVariableFVReal
    block = '${flow_blocks}'
  []
  #Porosity variable
  [porosity_var]
    type = MooseVariableFVReal
    block = '${flow_blocks}'
  []
  #Helium density variable
  [he_density]
    type = MooseVariableFVReal
  []
[]
#Applies variations to necessary AuxVariables
[AuxKernels]
  #Populates each block with the correct power density 
  [populate_power_source]
    type = FunctorElementalAux
    variable = 'power_source'
    functor = 'heat_source'
  []
  #Populates each block with the correct porosity
  [populate_porosity_var]
    type = FunctorElementalAux
    variable = 'porosity_var'
    functor = 'porosity'
  []
  #Computes the interstitial velocity from the superficial velocity 
  [populate_interstitial_vel]
    type = ParsedAux
    variable = 'interstitial_vel'
    coupled_variables = 'superficial_vel_z porosity_var'
    expression = 'superficial_vel_z / porosity_var'
  []
  #Computes helium density
  [compute_he_density]
    type = ADFunctorElementalAux
    variable = he_density
    functor = 'rho'
    execute_on = 'TIMESTEP_END'
  []
[]
#Gives relevant data about solution variables and physical parameters 
[Postprocessors]
  #Pressure drop in the HTTF
  [pdrop_total]
    type = PressureDrop
    pressure = pressure
    upstream_boundary = 'inlet_top'
    downstream_boundary = 'fuel_inner_bot'
    boundary = 'inlet_top fuel_inner_bot'
  []
  #Inlet pressure
  [p_in]
    type = FunctionSideAverage
    function = pressure_fn
    boundary = 'inlet_top'
  []
  [p_in_check]
    type = SideAverageValue
    variable = pressure
    boundary = 'inlet_top'
  []
  #Inlet mass flow rate
  [mdot_in]
    type = FunctionSideAverage
    function = mdot_fn
    boundary = 'inlet_top'
  []
  #Outlet coolant temperature
  [outlet_temp]
    type = SideAverageValue
    variable = T_fluid
    boundary = 'fuel_inner_bot fuel_outer_bot'
  []
  #Inlet coolant temperature
  [inlet_temp]
    type = SideAverageValue
    variable = T_fluid
    boundary = 'inlet_top'
  []
  #Temperature difference between the inlet and outlet boundaries
  [delta_T]
    type = ParsedPostprocessor
    pp_names = 'outlet_temp inlet_temp'
    function = 'outlet_temp - inlet_temp'
  []
  #Flow area 
  [flow_area]
    type = AreaPostprocessor
    boundary = 'inlet_top'
  []
  #Average velocity in the inlet
  [average_velocity]
    type = SideAverageValue
    boundary = 'inlet_top'
    variable = 'superficial_vel_z'
  []
  #Average density in the inlet
  [average_density]
    type = SideAverageValue
    boundary = 'inlet_top'
    variable = 'he_density'
  []
  #Area of the core
  [area_top]
    type = AreaPostprocessor
    boundary = 'inlet_top'
  []
  #Average mass flow rate
  [mass_flow_rate]
    type = ParsedPostprocessor
    pp_names = 'average_velocity average_density area_top'
    function = 'average_velocity * average_density * area_top'
  []
  #Total power in the core
  [total_power]
    type = ElementIntegralVariablePostprocessor
    variable = power_source
    block = 'core_iref_oring_cen core_iref_oring_cor core_fuel_inner 
    core_fuel_outer core_oref_iring_cor core_oref_iring_con'
  []
  [T_core_inner_top]
    type = ElementAverageValue
    variable = T_fluid
    block = '42 43'
  []
  [T_core_inner_middle]
    type = ElementAverageValue
    variable = T_fluid
    block = '12 13'
  []
  [T_core_inner_bottom]
    type = ElementAverageValue
    variable = T_fluid
    block = '32 33'
  []
  [T_core_middle_top]
    type = ElementAverageValue
    variable = T_fluid
    block = '44'
  []
  [T_core_middle_middle]
    type = ElementAverageValue
    variable = T_fluid
    block = '14'
  []
  [T_core_middle_bottom]
    type = ElementAverageValue
    variable = T_fluid
    block = '34'
  []
  [T_core_outer_top]
    type = ElementAverageValue
    variable = T_fluid
    block = '47'
  []
  [T_core_outer_middle]
    type = ElementAverageValue
    variable = T_fluid
    block = '17'
  []
  [T_core_outer_bottom]
    type = ElementAverageValue
    variable = T_fluid
    block = '37'
  []
  [T_core_outer_reflector]
    type = ElementAverageValue
    variable = T_solid
    block = '20'
  []
  [T_core_inner_reflector]
    type = ElementAverageValue
    variable = T_solid
    block = '11 41 31'
  []
  #Core fluid temperatures
  [TF-1504]
    type = PointValue
    variable =  T_fluid
    point = '${fparse R_l -0.026030002} 0.258540009 1.479'
    execute_on = 'TIMESTEP_END FINAL'
  [] 
  [TF-1520]
    type = PointValue
    variable = T_fluid
    point = '${fparse R_l + 0.338390015} 0.01503 1.479'
    execute_on = 'TIMESTEP_END FINAL'
  []
  [TS-1716]
    type = PointValue
    variable = T_solid
    point = '0.05206 -0.1297 1.876'
    execute_on = 'TIMESTEP_END FINAL'
  []
  [TS-1740]
    type = PointValue
    variable = T_solid
    point = '0.286329987 -0.33062 1.876'
    execute_on = 'TIMESTEP_END FINAL'
  []
  [TF-1708]
    type = PointValue
    variable = T_fluid
    point = '${fparse R_l + -0.286329987} 0.285540009 1.876'
    execute_on = 'TIMESTEP_END FINAL'
  []
  #Effective thermal conductivity imported from the subassembly 
  [k_eff]
    type = Receiver
    execute_on = 'TIMESTEP_END FINAL'
    default = 4.5
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -ksp_gmres_restart -pc_factor_shift_type '
    petsc_options_value = 'lu       20                  NONZERO              '
  []
[]


[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  start_time = 0
  end_time = 220000
  l_max_its = 100
  dtmax = 1000
  line_search = basic
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    iteration_window = 2
    optimal_iterations = 10
  []
  nl_abs_tol = 1e-1
  nl_rel_tol = 1e-4
  nl_max_its = 20
  nl_forced_its = 2
  automatic_scaling = true
  compute_scaling_once = false
  off_diagonals_in_auto_scaling = true 
  steady_state_detection = true
[]

[Outputs]
  exodus = true
  csv = true
  checkpoint = true
  print_linear_residuals = true 
[]

#================================================================================================
#                                          MultiApps
#================================================================================================
#Subscale models which communicate with main app porous medium model
[MultiApps]
  # Sub app for the subscale model of the homogenized region number 3 (no rotation)
  [subscale_3_1]
    type = TransientMultiApp
    positions = '0 0.18034 2.286'
    input_files = 'assembly_3_1.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
  # Sub app for the subscale model of the homogenized region number 3 (60 degrees rotation)
  [subscale_3_2]
    type = TransientMultiApp
    positions = '0.156179 0.09017 2.286'
    input_files = 'assembly_3_2.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 3 (120 degrees rotation)
  [subscale_3_3]
    type = TransientMultiApp
    positions = '0.156179 -0.09017 2.286'
    input_files = 'assembly_3_3.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 3 (180 degrees rotation)
  [subscale_3_4]
    type = TransientMultiApp
    positions = '0 -0.18034 2.286'
    input_files = 'assembly_3_4.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 3 (240 degrees rotation)
  [subscale_3_5]
    type = TransientMultiApp
    positions = '-0.156179 -0.09017 2.286'
    input_files = 'assembly_3_5.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 3 (300 degrees rotation)
  [subscale_3_6]
    type = TransientMultiApp
    positions = '-0.156179 0.09017 2.286'
    input_files = 'assembly_3_6.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 4
  [subscale_4]
    type = TransientMultiApp
    positions_file = 'positions_inner.txt'
    input_files = 'assembly_4.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 5 (60 degrees rotation)
  [subscale_5_1]
    type = TransientMultiApp
    positions_file = 'positions_5_1.txt'
    input_files = 'assembly_5_1.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 5 (120 degrees rotation)
  [subscale_5_2]
    type = TransientMultiApp
    positions_file = 'positions_5_2.txt'
    input_files = 'assembly_5_2.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 5 (180 degrees rotation)
  [subscale_5_3]
    type = TransientMultiApp
    positions_file = 'positions_5_3.txt'
    input_files = 'assembly_5_3.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
  # Sub app for the subscale model of the homogenized region number 5 (240 degrees rotation)
  [subscale_5_4]
    type = TransientMultiApp
    positions_file = 'positions_5_4.txt'
    input_files = 'assembly_5_4.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 5 (300 degrees rotation)
  [subscale_5_5]
    type = TransientMultiApp
    positions_file = 'positions_5_5.txt'
    input_files = 'assembly_5_5.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 5 (no rotation)
  [subscale_5_6]
    type = TransientMultiApp
    positions_file = 'positions_5_6.txt'
    input_files = 'assembly_5_6.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 6 (no rotation)
  [subscale_6_1]
    type = TransientMultiApp
    positions = '0 0.45085 2.286'
    input_files = 'assembly_6_1.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 6 (60 degrees rotation)
  [subscale_6_2]
    type = TransientMultiApp
    positions = '0.390448 0.225425 2.286'
    input_files = 'assembly_6_2.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 6 (120 degrees rotation)
  [subscale_6_3]
    type = TransientMultiApp
    positions = '0.390448 -0.225425 2.286'
    input_files = 'assembly_6_3.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 6 (180 degrees rotation)
  [subscale_6_4]
    type = TransientMultiApp
    positions = '0 -0.45085 2.286'
    input_files = 'assembly_6_4.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 6 (240 degrees rotation)
  [subscale_6_5]
    type = TransientMultiApp
    positions = '-0.390448 -0.225425 2.286'
    input_files = 'assembly_6_5.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 6 (300 degrees rotation)
  [subscale_6_6]
    type = TransientMultiApp
    positions = '-0.390448 0.225425 2.286'
    input_files = 'assembly_6_6.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 7 (no rotation)
  [subscale_7_1]
    type = TransientMultiApp
    positions_file = 'positions_7_1.txt'
    input_files = 'assembly_7_1.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 7 (60 degrees rotation)
  [subscale_7_2]
    type = TransientMultiApp
    positions_file = 'positions_7_2.txt'
    input_files = 'assembly_7_2.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 7 (120 degrees rotation)
  [subscale_7_3]
    type = TransientMultiApp
    positions_file = 'positions_7_3.txt'
    input_files = 'assembly_7_3.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 7 (180 degrees rotation)
  [subscale_7_4]
    type = TransientMultiApp
    positions_file = 'positions_7_4.txt'
    input_files = 'assembly_7_4.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 7 (240 degrees rotation)
  [subscale_7_5]
    type = TransientMultiApp
    positions_file = 'positions_7_5.txt'
    input_files = 'assembly_7_5.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
   # Sub app for the subscale model of the homogenized region number 7 (300 degrees rotation)
  [subscale_7_6]
    type = TransientMultiApp
    positions_file = 'positions_7_6.txt'
    input_files = 'assembly_7_6.i' 
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    execute_on = ' TIMESTEP_END'
  []
[]
#Transfers between porous medium model and subscale models 
[Transfers]
  #transfers the effective thermal conductivity from the SubApp to the mainApp
  [k_eff_transfer]
    type = MultiAppPostprocessorTransfer
    from_multi_app = subscale_4
    from_postprocessor = k_eff
    to_postprocessor = k_eff
    reduction_type = average 
  []
  #Transfer of pressure from the MainApp to the subApp
  [pressure_to_subscale_3_1]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_3_1
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_3_2]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_3_2
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_3_3]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_3_3
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_3_4]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_3_4
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_3_5]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_3_5
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_3_6]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_3_6
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_4]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_4
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_5_1]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_5_1
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_5_2]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_5_2
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_5_3]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_5_3
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_5_4]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_5_4
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_5_5]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_5_5
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_5_6]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_5_6
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_6_1]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_6_1
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_6_2]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_6_2
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_6_3]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_6_3
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_6_4]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_6_4
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_6_5]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_6_5
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_6_6]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_6_6
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_7_1]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_7_1
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_7_2]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_7_2
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_7_3]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_7_3
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_7_4]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_7_4
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_7_5]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_7_5
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  [pressure_to_subscale_7_6]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_7_6
    from_postprocessor = p_in
    to_postprocessor = system_pressure
  []
  #Transfer of the mass flow rate from the mainApp to the subApps
  [mdot_to_subscale_3_1]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_3_1
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_3_2]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_3_2
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_3_3]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_3_3
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_3_4]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_3_4
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_3_5]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_3_5
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_3_6]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_3_6
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_4]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_4
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_5_1]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_5_1
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_5_2]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_5_2
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_5_3]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_5_3
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_5_4]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_5_4
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_5_5]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_5_5
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_5_6]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_5_6
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_6_1]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_6_1
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_6_2]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_6_2
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_6_3]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_6_3
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_6_4]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_6_4
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_6_5]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_6_5
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_6_6]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_6_6
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_7_1]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_7_1
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_7_2]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_7_2
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_7_3]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_7_3
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_7_4]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_7_4
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_7_5]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_7_5
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  [mdot_to_subscale_7_6]
    type = MultiAppPostprocessorTransfer
    to_multi_app = subscale_7_6
    from_postprocessor = mdot_in
    to_postprocessor = m_dot_total
  []
  #Transfer of the solid temperatures to the subscale model on the boundaries
  [T_boundary_to_subscale_3_1]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_3_1
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '171 172 173 174 175 176' # Boundaries of the assembly (not including the gaps and the heater rods)
  []
  [T_boundary_to_subscale_3_2]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_3_2
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '171 172 173 174 175 176'
  []
  [T_boundary_to_subscale_3_3]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_3_3
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '171 172 173 174 175 176'
  []
  [T_boundary_to_subscale_3_4]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_3_4
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '171 172 173 174 175 176'
  []
  [T_boundary_to_subscale_3_5]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_3_5
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '171 172 173 174 175 176'
  []
  [T_boundary_to_subscale_3_6]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_3_6
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '171 172 173 174 175 176'
  []
  [T_boundary_to_subscale_4]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_4
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '161 162 163 164 165 166' # Boundaries of the assembly (not including the gaps and the heater rods)
  []
  [T_boundary_to_subscale_5_1]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_5_1
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '151 152 153 154 155 156' # Boundaries of the assembly (not including the gaps and the heater rods)
  []
  [T_boundary_to_subscale_5_2]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_5_2
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '151 152 153 154 155 156'
  []
  [T_boundary_to_subscale_5_3]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_5_3
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '151 152 153 154 155 156'
  []
  [T_boundary_to_subscale_5_4]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_5_4
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '151 152 153 154 155 156'
  []
  [T_boundary_to_subscale_5_5]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_5_5
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '151 152 153 154 155 156'
  []
  [T_boundary_to_subscale_5_6]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_5_6
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '151 152 153 154 155 156'
  []
  [T_boundary_to_subscale_6_1]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_6_1
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '141 142 143 144 145 146' # Boundaries of the assembly (not including the gaps and the heater rods)
  []
  [T_boundary_to_subscale_6_2]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_6_2
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '141 142 143 144 145 146'
  []
  [T_boundary_to_subscale_6_3]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_6_3
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '141 142 143 144 145 146'
  []
  [T_boundary_to_subscale_6_4]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_6_4
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '141 142 143 144 145 146'
  []
  [T_boundary_to_subscale_6_5]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_6_5
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '141 142 143 144 145 146'
  []
  [T_boundary_to_subscale_6_6]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_6_6
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '141 142 143 144 145 146'
  []
  [T_boundary_to_subscale_7_1]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_7_1
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '181 182 183 184 185 186' # Boundaries of the assembly (not including the gaps and the heater rods)
  []
  [T_boundary_to_subscale_7_2]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_7_2
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '181 182 183 184 185 186'
  []
  [T_boundary_to_subscale_7_3]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_7_3
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '181 182 183 184 185 186'
  []
  [T_boundary_to_subscale_7_4]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_7_4
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '181 182 183 184 185 186'
  []
  [T_boundary_to_subscale_7_5]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_7_5
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '181 182 183 184 185 186'
  []
  [T_boundary_to_subscale_7_6]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = subscale_7_6
    source_variable = T_solid
    variable = T_boundary
    to_boundaries = '181 182 183 184 185 186'
  []
[]

[Debug]
  show_var_residual_norms = true
  check_jacobian = true
[]
