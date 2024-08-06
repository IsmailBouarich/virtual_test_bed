# ==============================================================================
# High Temperature Transient Facility (HTTF)
# Subscale solve of the HTTF steady state for assembly number 4 in the hexagonal homogenzation scheme
# ------------------------------------------------------------------------------
# Idaho Falls, INL, 08/2024
# Author(s): Ismail Bouarich, Dr. Lise Charlot
# If using or referring to this model, please cite as explained in
# https://mooseframework.inl.gov/virtual_test_bed/citing.html
# ==============================================================================

!include common.i
# Description of hexagonal assembly number 4
radius_large_coolant = ${units 7.9375 mm -> m}

assembly_apothem = ${units 75.4 mm -> m}
pincell_apothem = ${units 15.08 mm -> m}
assembly_height = 1.9812 #m
T_in = 500 #K

fuel_density = 1750 #kg/m3
ns = 4
#================================================================================================
#                                          Meshing
#================================================================================================
#Creates mesh for assembly 4
[Mesh]

  #Pincells
  #Fuel rod pincells
  [fuel_rod_pincell_1]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 2.
    background_block_ids = '101'
    polygon_size = ${pincell_apothem}
    polygon_size_style = 'apothem'
    ring_radii = '0.0076 ${radius_fuel}'
    ring_intervals = '2 1'
    ring_block_ids = '1 111 11'
    preserve_volumes = on
    quad_center_elements = false
  []
  [fuel_rod_pincell_2]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 2.
    background_block_ids = '102'
    polygon_size = ${pincell_apothem}
    polygon_size_style = 'apothem'
    ring_radii = '0.0076 ${radius_fuel}'
    ring_intervals = '1 1'
    ring_block_ids = '2 12'
    preserve_volumes = on
    quad_center_elements = false
  []
  [fuel_rod_pincell_3]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 2.
    background_block_ids = '103'
    polygon_size = ${pincell_apothem}
    polygon_size_style = 'apothem'
    ring_radii = '0.0076 ${radius_fuel}'
    ring_intervals = '1 1'
    ring_block_ids = '3 13'
    preserve_volumes = on
    quad_center_elements = false
  []
  [fuel_rod_pincell_4]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 2.
    background_block_ids = '104'
    polygon_size = ${pincell_apothem}
    polygon_size_style = 'apothem'
    ring_radii = '0.0076 ${radius_fuel}'
    ring_intervals = '1 1'
    ring_block_ids = '4 14'
    preserve_volumes = on
    quad_center_elements = false
  []
  [fuel_rod_pincell_5]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 2.
    background_block_ids = '105'
    polygon_size = ${pincell_apothem}
    polygon_size_style = 'apothem'
    ring_radii = '0.0076 ${radius_fuel}'
    ring_intervals = '1 1'
    ring_block_ids = '5 15'
    preserve_volumes = on
    quad_center_elements = false
  []
  [fuel_rod_pincell_6]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 2.
    background_block_ids = '106'
    polygon_size = ${pincell_apothem}
    polygon_size_style = 'apothem'
    ring_radii = '0.0076 ${radius_fuel}'
    ring_intervals = '1 1'
    ring_block_ids = '6 16'
    preserve_volumes = on
    quad_center_elements = false
  []
  [fuel_rod_pincell_7]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 2.
    background_block_ids = '107'
    polygon_size = ${pincell_apothem}
    polygon_size_style = 'apothem'
    ring_radii = '0.0076 ${radius_fuel}'
    ring_intervals = '1 1'
    ring_block_ids = '7 17'
    preserve_volumes = on
    quad_center_elements = false
  []
  #Large coolant pincell
  [large_coolant_pincell]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6 # must be six to use hex pattern
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 2.
    background_block_ids = '108'
    polygon_size = ${pincell_apothem}
    polygon_size_style = 'apothem'
    ring_radii = '${radius_large_coolant}'
    ring_intervals ='1'
    ring_block_ids = '8'
    ring_block_names = 'large_coolant_pincell'
    preserve_volumes = on
    quad_center_elements = false
  []
  #Dummy pincell used for hexagonal assembly that will be deleted afterwards
  [dummy_pincell]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6 # must be six to use hex pattern
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 2.
    background_block_ids = '109'
    polygon_size = ${pincell_apothem}
    polygon_size_style = 'apothem'
    ring_radii = '${radius_large_coolant}'
    ring_intervals ='1'
    ring_block_ids = '9'
    ring_block_names = 'dummy_pincell'
    preserve_volumes = on
    quad_center_elements = false
  []

  ### Assembly

  [assembly]
    type = PatternedHexMeshGenerator
    inputs = 'fuel_rod_pincell_1 fuel_rod_pincell_2 fuel_rod_pincell_3 fuel_rod_pincell_4 fuel_rod_pincell_5 fuel_rod_pincell_6 fuel_rod_pincell_7 large_coolant_pincell dummy_pincell'
    # Pattern ID     0                 1                   2                  3                   4                 5                   6                   7                 8
    hexagon_size = ${assembly_apothem}
    pattern = '  8 2 8;
                1 7 7 3;
               8 7 0 7 8;
                6 7 7 4;
                 8 5 8;'
  []
  #Peripheral trimming 
  [trim]
    type = HexagonMeshTrimmer
    input = assembly
    trim_peripheral_region = '1 1 1 1 1 1'
    peripheral_trimming_section_boundary = peripheral_section
  []
  #Assembly rotation to match the HTTF core assemblies
  [rotate_assembly]
    type = TransformGenerator
    input = trim
    transform = ROTATE
    vector_value = '0 0 90'
  []
  #Extrusion to 3D
  [extrude]
    type = AdvancedExtruderGenerator
    input = rotate_assembly
    heights = '${assembly_height}'
    num_layers = '50'
    direction = '0 0 -1'
  []
  #Deletion of dummy pincells
  [dummy_delete]
    type = BlockDeletionGenerator
    block = '9 109'
    input = extrude
  []
  #Deletion of unwanted subdomain to match the desired assembly
  [plane_deletion_1]
    type = PlaneDeletionGenerator
    input = dummy_delete
    point = '0.04524 0 0'
    normal = '1.0 0 0'
  []

  [plane_deletion_2]
    type = PlaneDeletionGenerator
    input = plane_deletion_1
    point = '0 0.05224 0'
    normal = '0.02612 0.04524 0'
  []

  [plane_deletion_3]
    type = PlaneDeletionGenerator
    input = plane_deletion_2
    point = '0 0.05224 0'
    normal = '-0.02612 0.04524 0'
  []

  [plane_deletion_4]
    type = PlaneDeletionGenerator
    input = plane_deletion_3
    point = '-0.04524 0 0'
    normal = '-1 0 0'
  []

  [plane_deletion_5]
    type = PlaneDeletionGenerator
    input = plane_deletion_4
    point = '0 -0.05224 0'
    normal = '-0.02612 -0.04524 0'
  []

  [plane_deletion_6]
    type = PlaneDeletionGenerator
    input = plane_deletion_5
    point = '0 -0.05224 0'
    normal = '0.02612 -0.04524 0'
  []
  #Adds coolant boudnary
  [add_coolant_boundary]
    type = SideSetsBetweenSubdomainsGenerator
    input = plane_deletion_6
    primary_block = '108'
    paired_block = '8'
    new_boundary = coolant_boundary
  []
  #Deletion of the solid block to create coolant gaps
  [assembly_final]
    type = BlockDeletionGenerator
    block = '8'
    input = add_coolant_boundary
  []
  #Adds assembly boundaries without including the heater rods and their gaps 
  [assembly_boundary]
    type = SideSetsFromNormalsGenerator
    input = assembly_final
    new_boundary = '161 162 163 164 165 166'
    normals = '1 0 0 0.02612 0.04524 0 -0.02612 0.04524 0 -1 0 0 -0.02612 -0.04524 0 0.02612 -0.04524 0'
    included_subdomains = '101 102 103 104 105 106 107 108'
  []
  #Definition of outer gaps 
  [gap_outer_1]
    type = SideSetsBetweenSubdomainsGenerator
    input = assembly_boundary
    primary_block = '101'
    paired_block = '11'
    new_boundary = gap_outer_1
  []
  [gap_outer_2]
    type = SideSetsBetweenSubdomainsGenerator
    input = gap_outer_1
    primary_block = '102'
    paired_block = '12'
    new_boundary = gap_outer_2
  []
  [gap_outer_3]
    type = SideSetsBetweenSubdomainsGenerator
    input = gap_outer_2
    primary_block = '103'
    paired_block = '13'
    new_boundary = gap_outer_3
  []
  [gap_outer_4]
    type = SideSetsBetweenSubdomainsGenerator
    input = gap_outer_3
    primary_block = '104'
    paired_block = '14'
    new_boundary = gap_outer_4
  []
  [gap_outer_5]
    type = SideSetsBetweenSubdomainsGenerator
    input = gap_outer_4
    primary_block = '105'
    paired_block = '15'
    new_boundary = gap_outer_5
  []
  [gap_outer_6]
    type = SideSetsBetweenSubdomainsGenerator
    input = gap_outer_5
    primary_block = '106'
    paired_block = '16'
    new_boundary = gap_outer_6
  []
  [gap_outer_7]
    type = SideSetsBetweenSubdomainsGenerator
    input = gap_outer_6
    primary_block = '107'
    paired_block = '17'
    new_boundary = gap_outer_7
  []
  #Definition of inner gaps 
  [gap_inner_1]
    type = SideSetsBetweenSubdomainsGenerator
    input = gap_outer_7
    primary_block = '111'
    paired_block = '11'
    new_boundary = gap_inner_1
  []
  [gap_inner_2]
    type = SideSetsBetweenSubdomainsGenerator
    input = gap_inner_1
    primary_block = '2'
    paired_block = '12'
    new_boundary = gap_inner_2
  []
  [gap_inner_3]
    type = SideSetsBetweenSubdomainsGenerator
    input = gap_inner_2
    primary_block = '3'
    paired_block = '13'
    new_boundary = gap_inner_3
  []
  [gap_inner_4]
    type = SideSetsBetweenSubdomainsGenerator
    input = gap_inner_3
    primary_block = '4'
    paired_block = '14'
    new_boundary = gap_inner_4
  []
  [gap_inner_5]
    type = SideSetsBetweenSubdomainsGenerator
    input = gap_inner_4
    primary_block = '5'
    paired_block = '15'
    new_boundary = gap_inner_5
  []
  [gap_inner_6]
    type = SideSetsBetweenSubdomainsGenerator
    input = gap_inner_5
    primary_block = '6'
    paired_block = '16'
    new_boundary = gap_inner_6
  []
  [gap_inner_7]
    type = SideSetsBetweenSubdomainsGenerator
    input = gap_inner_6
    primary_block = '7'
    paired_block = '17'
    new_boundary = gap_inner_7
  []
  [rotate_assembly_2]
    type = TransformGenerator
    input = gap_inner_7
    transform = ROTATE
    vector_value = '0 0 30'
  []
  #Deletion of the solid block to create heater gaps
  [assembly_final_2]
    type = BlockDeletionGenerator
    block = '11 12 13 14 15 16 17'
    input = rotate_assembly_2
  []
[]

#================================================================================================
#                                          Heat Conduction
#================================================================================================
#Creates temperature variable
[Variables]
  [T]
    initial_condition = 500
  []
[]

[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
  # allow_initial_conditions_with_restart = true
  restart_file_base = HTTF_steady_state_out_subscale_426_cp/LATEST
[]
#Creates relevant AuxVariables such as fluid temperature, power density...
[AuxVariables]
  #Fluid temperature
  [T_fluid]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 500
  []
  [power_density]
    family = MONOMIAL
    order = CONSTANT
    block = '1 2 3 4 5 6 7 111'
  []
  [fuel_density]
  []
  #Heat transfer coefficient
  [htc]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 200
  []
  #Wall temperature of the coolant gaps
  [T_wall]
    initial_condition = 500
  []
  #Boundary temperature in the assembly boundaries
  [T_boundary]
    initial_condition = 500
  []
[]
#Dictates AuxVariables evolution
[AuxKernels]
  [density_fuel_ring]
    type = ConstantAux
    variable = fuel_density
    block = '1 2 3 4 5 6 7 111'
  []
  #Populates with a spatial value returned from wall temperature UserObject 
  [T_wall_aux]
    type = SpatialUserObjectAux
    variable = T_wall
    user_object = T_wall_coolant_uo
  []
  #Populates the power density following a specified function
  [power_density_aux]
    type = FunctionAux
    function = power_density_sub
    variable = power_density
    block = '1 2 3 4 5 6 7 111'
  []
[]
#Heat conduction equation modeling
[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = T
  []
  [time_derivative]
    type = HeatConductionTimeDerivative
    variable = T
  []
  #Applies the power as heat source in the heater rods
  [heat_source_fuel]
    type = CoupledForce
    variable = T
    block = '1 2 3 4 5 6 7 111'
    v = power_density
  []
[]
#Gap heat transfer modeling 
[ThermalContact]
  #Gap heat transfer between each inner gap and its corresponding outer gap
  [gap_heat_transfer_1]
    type = GapHeatTransfer
    primary = 'gap_inner_1'
    secondary = 'gap_outer_1'
    quadrature = true
    gap_geometry_type = CYLINDER
    cylinder_axis_point_1 = '0 0 ${fparse -assembly_height}'
    cylinder_axis_point_2 = '0 0 0'
    emissivity_primary = 0.9
    emissivity_secondary = 0.581
    variable = T
    gap_conductivity = 0.16
  []
  [gap_heat_transfer_2]
    type = GapHeatTransfer
    primary = 'gap_inner_2'
    secondary = 'gap_outer_2'
    quadrature = true
    gap_geometry_type = CYLINDER
    cylinder_axis_point_1 = '0.0522387 0 ${fparse -assembly_height}'
    cylinder_axis_point_2 = '0.0522387 0 0'
    emissivity_primary = 0.9
    emissivity_secondary = 0.581
    variable = T
    gap_conductivity = 0.16
  []
  [gap_heat_transfer_3]
    type = GapHeatTransfer
    primary = 'gap_inner_3'
    secondary = 'gap_outer_3'
    quadrature = true
    gap_geometry_type = CYLINDER
    cylinder_axis_point_1 = '0.0261213 -0.0452388 ${fparse -assembly_height}'
    cylinder_axis_point_2 = '0.0261213 -0.0452388 0'
    emissivity_primary = 0.9
    emissivity_secondary = 0.581
    variable = T
    gap_conductivity = 0.16
  []
  [gap_heat_transfer_4]
    type = GapHeatTransfer
    primary = 'gap_inner_4'
    secondary = 'gap_outer_4'
    quadrature = true
    gap_geometry_type = CYLINDER
    cylinder_axis_point_1 = '-0.0261213 -0.0452388 ${fparse -assembly_height}'
    cylinder_axis_point_2 = '-0.0261213 -0.0452388 0'
    emissivity_primary = 0.9
    emissivity_secondary = 0.581
    variable = T
    gap_conductivity = 0.16
  []
  [gap_heat_transfer_5]
    type = GapHeatTransfer
    primary = 'gap_inner_5'
    secondary = 'gap_outer_5'
    quadrature = true
    gap_geometry_type = CYLINDER
    cylinder_axis_point_1 = '-0.0522387 0 ${fparse -assembly_height}'
    cylinder_axis_point_2 = '-0.0522387 0 0'
    emissivity_primary = 0.9
    emissivity_secondary = 0.581
    variable = T
    gap_conductivity = 0.16
  []
  [gap_heat_transfer_6]
    type = GapHeatTransfer
    primary = 'gap_inner_6'
    secondary = 'gap_outer_6'
    quadrature = true
    gap_geometry_type = CYLINDER
    cylinder_axis_point_1 = '-0.0261213 0.0452388 ${fparse -assembly_height}'
    cylinder_axis_point_2 = '-0.0261213 0.0452388 0'
    emissivity_primary = 0.9
    emissivity_secondary = 0.581
    variable = T
    gap_conductivity = 0.16
  []
  [gap_heat_transfer_7]
    type = GapHeatTransfer
    primary = 'gap_inner_7'
    secondary = 'gap_outer_7'
    quadrature = true
    gap_geometry_type = CYLINDER
    cylinder_axis_point_1 = '0.0261213 0.0452388 ${fparse -assembly_height}'
    cylinder_axis_point_2 = '0.0261213 0.0452388 0'
    emissivity_primary = 0.9
    emissivity_secondary = 0.581
    variable = T
    gap_conductivity = 0.16
  []
[]

#Implements material properties in the system
[Materials]
  [ceramics_thermal]
    type = HeatConductionMaterial
    block = '101 102 103 104 105 106 107 108'
    temp = T
    thermal_conductivity_temperature_function = k_greencast_fn
    specific_heat_temperature_function = cp_greencast_fn
  []
  [ceramics_density]
    type = Density
    block = '101 102 103 104 105 106 107 108'
    density = '2912'
  []
  [fuel_thermal]
    type = HeatConductionMaterial
    block = '1 2 3 4 5 6 7 111'
    temp = T
    thermal_conductivity_temperature_function = k_G348
    specific_heat_temperature_function = cp_G348
  []
  [fuel_density]
    type = Density
    block = '1 2 3 4 5 6 7 111'
    density = '1860'
  []
[]
#Assembly boundary conditions
[BCs]
  # Couples far field fluid temperature to the temperature of the solid with a heater transfer coefficient htc
  [cooling_channels]
    type = CoupledConvectiveHeatFluxBC
    boundary = coolant_boundary
    T_infinity = T_fluid
    htc = htc
    variable = T
  []
  # Matches the value of the temperature T in assembly boundaries to the value of the solid temperature in the porous media model
  [boundary_heat]
    type = MatchedValueBC
    boundary = '161 162 163 164 165 166'
    v = T_boundary
    # htc = 1e5
    variable = T
  []
[]

#================================================================================================
#                                          Post Processors
#================================================================================================
# UserObjects used in transfer to THM subscale 
[UserObjects]
  # Determines wall temperature for the coolant boundary and applies them to the positions listed in 'large_coolants_positions.txt'
  [T_wall_coolant_uo]
    type = NearestPointLayeredSideAverage
    variable = T
    boundary = coolant_boundary
    direction = z
    points_file = large_coolants_positions.txt
    num_layers = 50
  []
[]
#Gives relevant data about solution variables and physical parameters 
[Postprocessors]
  #Assembly boundary average temperature 
  [temp_boundary_subscale]
    type = SideAverageValue
    variable = T_boundary
    boundary = '161 162 163 164 165 166'
  []
  #Coolant channel average wall temperature
  [temp_wall_channel]
    type = SideAverageValue
    variable = T_wall
    boundary = 'coolant_boundary'
  []
  #Middle inner gap average temperature
  [temp_gap_inner_middle]
    type = SideAverageValue
    variable = T
    boundary = 'gap_inner_1'
  []
  #Middle outer gap average temperature
  [temp_gap_outer_middle]
    type = SideAverageValue
    variable = T
    boundary = 'gap_outer_1'
  []
  #Average heat transfer coefficient value
  [ht_value]
    type = ElementAverageValue
    variable = htc
  []
  #Fluid temperature
  [T_fluid_subscale]
    type = ElementAverageValue
    variable = T_fluid
  []
  [total_power_subscale]
    type = ElementIntegralVariablePostprocessor
    variable = power_density
    block = '1 2 3 4 5 6 7 111'
  []
  [m_dot_total]
    type = Receiver
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    default = 1
  []
  [system_pressure]
    type = Receiver
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    default = 7e5
  []
  #Amount of heat transferred to coolant channels
  [heat_to_channels]
    type = ConvectiveHeatTransferSideIntegral
    T_solid = T
    T_fluid_var = T_fluid
    htc_var = htc
    boundary = coolant_boundary
  []
  #Heat flux in coolant boundary
  [heat_flux]
    type = SideFluxIntegral
    variable = T
    boundary = '161 162 163 164 165 166'
    diffusivity = thermal_conductivity
  []
  #Effective thermal conductivity
  [k_eff]
    type = ParsedPostprocessor
    expression = 'heat_flux * (0.04524/(temp_gap_outer_middle-temp_boundary_subscale))'
    pp_names = 'heat_flux temp_gap_outer_middle temp_boundary_subscale'
  []
[]

#================================================================================================
#                                          Execution
#================================================================================================

[Executioner]
  type = Transient
  start_time = 0
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
  []
  dtmax = 150000
  dtmin =1e-4
  end_time = 2200000

  line_search = basic
  solve_type = NEWTON

  automatic_scaling = true

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  fixed_point_abs_tol = 1e-2
  fixed_point_max_its = 25
  fixed_point_rel_tol = 1e-4
  nl_rel_tol = 1e-2
  nl_abs_tol = 1e-4
  nl_max_its = 40
[]

[Outputs]
  # [out]
  #   type = Exodus
  # []
  csv = true 
  # checkpoint = true
[]

#================================================================================================
#                                          MultiApps
#================================================================================================
#THM (Thermal Hydraulics) channels subapp representing coolant flows
[MultiApps]
  [coolant_flow]
    type = TransientMultiApp
    positions_file = 'large_coolants_positions.txt'
    input_files = 'large_coolant_channel.i'
    max_procs_per_app = 1
    output_in_position = true
    sub_cycling = true
    output_sub_cycles = true
    max_catch_up_steps = 1e4
    max_failures = 10000
    # execute_on = ' TIMESTEP_END'
    # bounding_box_padding = '0.01 0.01 0.01'
    # ignore_solve_not_converge = true
    # keep_solution_during_restore = true
  []
[]
#Transfers between main app assembly and coolant channels subapp
[Transfers]
  #Transfer of wall temperature to subApp 1D coolant channel
  [T_wall_to_coolant]
    type = MultiAppGeneralFieldUserObjectTransfer
    source_user_object = T_wall_coolant_uo
    variable = T_wall_channel
    to_multi_app = coolant_flow
  []
  #Transfer of coolant temperature from the subApp (coolant channel) to the mainApp 
  [T_coolant]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = coolant_flow
    source_variable = T
    variable = T_fluid
  []
  #Transfer of heat transfer coefficient from the subApp (coolant channel) to the mainApp 
  [htc_from_coolant]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = coolant_flow
    source_variable = htc_coolant
    variable = htc
  []
  #Transfer of mass flow rate to subApp 1D coolant channel
  [mdot_transfer]
    type = MultiAppPostprocessorTransfer
    to_multi_app = coolant_flow
    from_postprocessor = m_dot_total
    to_postprocessor = m_dot_total
  []
  #Transfer of pressure to subApp 1D coolant channel
  [pressure_transfer]
    type = MultiAppPostprocessorTransfer
    to_multi_app = coolant_flow
    from_postprocessor = system_pressure
    to_postprocessor = system_pressure
  []
[]