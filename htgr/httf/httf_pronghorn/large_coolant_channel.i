# Description of coolant channels
radius_large_coolant = ${units 7.9375 mm -> m}
A_channel = '${fparse pi * radius_large_coolant * radius_large_coolant}'
P_hf_channel = '${fparse  2 * pi * radius_large_coolant}'
T_in = 650 # K
p_out = 7e5 # Pa
n_large_channels = 324
n_small_channels = 96
n_medium_channels = 96
n_large_bypass = 6
n_small_bypass = 36
D_l = 0.016
A_l = ${fparse ((D_l*D_l)/4)*pi*n_large_channels}
D_m = 0.013
A_m = ${fparse ((D_m*D_m)/4)*pi*n_medium_channels}
D_s = 0.01
A_s = ${fparse ((D_s*D_s)/4)*pi*n_small_channels}
D_lb = 0.019
A_lb = ${fparse ((D_lb*D_lb)/4)*pi*n_large_bypass}
D_sb = 0.016
A_sb = ${fparse ((D_sb*D_sb)/4)*pi*n_small_bypass}
A_total = ${fparse A_l + A_m + A_s + A_lb + A_sb}
m_dot_total = 1
m_dot_in = ${fparse (m_dot_total*(A_l/A_total))/n_large_channels}
R_He = 2076.9 # J/kg-K
rho_He = ${fparse p_out/(R_He*T_in)}
vel_in = ${fparse m_dot_in/(rho_He*A_channel)}
l_core = 2.864 # m
assembly_height = 1.9812 #m
[GlobalParams]
    initial_vel = ${vel_in} #m/s
    initial_T = ${T_in}
    initial_p = ${p_out}
    rdg_slope_reconstruction = full
    scaling_factor_1phase = '1 1e-2 1e-5'
    fp = helium
    closures = thm
[]
[Problem]
    allow_initial_conditions_with_restart = true
[]
[FluidProperties]
    [helium]
        type = IdealGasFluidProperties
        molar_mass = 4e-3
        gamma = 1.67
        k = 0.2556
        mu = 3.22639e-5
    []
[]
[Closures]
  [thm]
    type = Closures1PhaseTHM
  []
[]
[AuxVariables]
    [htc_coolant]
      initial_condition = 200
      family = MONOMIAL
      order = CONSTANT
    []
    [T_wall_channel]
        initial_condition = ${T_in}
    []
[]
[AuxKernels]
    [htc_aux]
        type = ADMaterialRealAux
        property = Hw
        variable = htc_coolant
    []
[]
[Components]
    [inlet_coolant_chan]
        type = InletMassFlowRateTemperature1Phase
        m_dot = ${m_dot_in}
        T = ${T_in}
        input = coolant_chan:in
    []
    [coolant_chan]
        type = FlowChannel1Phase
        position = '0 0 0'
        orientation = '0 0 -1'
        length = ${assembly_height}
        n_elems = 50
        A = ${A_channel}
        D_h = '${fparse 2 * radius_large_coolant}'
    []
    [outlet_chan]
        type = Outlet1Phase
        p = ${p_out}
        input = coolant_chan:out
    []
    [ht]
        type = HeatTransferFromExternalAppTemperature1Phase
        flow_channel = coolant_chan
        P_hf = ${P_hf_channel}
        T_ext = T_wall_channel
    []
[]
[ControlLogic]
  [mdot_control]
    type = ParsedFunctionControl
    function = 'if(t>0, m_dot_in, ${m_dot_in})'
    vals = 'm_dot_in'
    vars = 'm_dot_in'
  []
  [mdot_applied]
    type = SetComponentRealValueControl
    component = inlet_coolant_chan
    parameter = m_dot
    value = mdot_control:value
  []
  [inlet_temp_logic]
    type = ParsedFunctionControl
    function = 'if(t > 0, inlet_temp, ${T_in})'
    vals = 'inlet_temp'
    vars = 'inlet_temp'
  []
  [inlet_temp_crtl]
    type = SetComponentRealValueControl
    component = inlet_coolant_chan
    parameter = T
    value = inlet_temp_logic:value
  []
  [sys_pressure_logic]
    type = ParsedFunctionControl
    function = 'if(t > 0, system_pressure, ${p_out})'
    vals = 'system_pressure'
    vars = 'system_pressure'
  []
  [sys_pressure_ctrl]
    type = SetComponentRealValueControl
    component = outlet_chan
    parameter = p
    value = sys_pressure_logic:value
  []
[]
#================================================================================================
#                                          Post Processors
#================================================================================================
[Postprocessors]
    [heat_to_coolants]
        type = ADHeatRateConvection1Phase
        P_hf = ${P_hf_channel}
        T_wall = T_wall_channel
    []
    [m_dot_in_check]
        type = ADFlowBoundaryFlux1Phase
        boundary = 'inlet_coolant_chan'
        equation = mass
    []
    [m_dot_total]
        type = Receiver
        default = ${m_dot_total}
       execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [m_dot_in]
        type = ParsedPostprocessor
        expression = 'm_dot_total*(${A_l}/${A_total})/${n_large_channels}'
        pp_names = 'm_dot_total'
    []
    [T_out]
        type = SideAverageValue
        variable = T
        boundary = 'outlet_chan'
    []
    [p_in]
        type = SideAverageValue
        variable = p
        boundary = inlet_coolant_chan
    []
    [p_out]
        type = SideAverageValue
        variable = p
        boundary = outlet_chan
    []
    [delta_p]
        type = ParsedPostprocessor
        pp_names = 'p_in p_out'
        function = 'p_in - p_out'
    []
    [Hw_value]
        type = ADElementAverageMaterialProperty
        mat_prop = Hw
    []
    [T_wall_channel]
        type = ElementAverageValue
        variable = T_wall_channel
    []
    [inlet_temp]
      type = Receiver
      execute_on = 'INITIAL TIMESTEP_BEGIN'
      default = ${T_in}
    []
    [system_pressure]
      type = Receiver
      execute_on = 'INITIAL TIMESTEP_BEGIN'
      default = ${p_out}
    []
[]
[Executioner]
    type = Transient
    # start_time = 0
    [TimeStepper]
        type = IterationAdaptiveDT
        dt = 0.01
    []
    dtmin = 1e-4
    dtmax = 15000
    end_time = 220000

    solve_type = NEWTON
    petsc_options_iname = '-pc_type'
    petsc_options_value = 'lu'
    line_search = basic
    # fixed_point_abs_tol = 1e-2
    # fixed_point_max_its = 25
    # fixed_point_rel_tol = 1e-4
    nl_rel_tol = 1e-2
    nl_abs_tol = 1e-4
    nl_max_its = 40
[]
[Outputs]
    exodus = true
[]