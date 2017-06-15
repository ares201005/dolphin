[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 0.003
  ymin = 0
  ymax = 1.0
  nx = 10
  ny = 100
  elem_type = QUAD9
[]


[Variables]
  [./eqpot]
    order = FIRST
    family = LAGRANGE
    scaling = 1.e-2
  [../]

  [./potential]
    order = FIRST
    family = LAGRANGE
  [../]

  [./conA]
    order = FIRST
    family = LAGRANGE
    scaling = 1.e-1
  [../]

  [./conC]
    order = FIRST
    family = LAGRANGE
    scaling = 1.e-1
  [../]
[]

# The Preconditioning block
[Preconditioning]
  active = 'mySMP'
  [./mySMP]
    type = SMP
    full = true 
  [../]
[]

[AuxVariables]
  [./totConA]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./totConC]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./eqpot_diff]
    type = Diffusion
    variable = eqpot
  [../]
  [./eqpot_induced]
    type = EqPotential
    variable = eqpot
    dielectric = 138801
    coefficient = 38.69565  # z_i  e  /KT
    bulkconc = 1
  [../]

  [./potential_diff]
    type = Diffusion
    variable = potential
  [../]

  [./potential_induced_A]
     type = Potential
     variable = potential
     concentration = conA
     dielectric = -138801  #0.1388015 # inverse of dielectric constant, z_i  e  /epsilon
  [../]

  [./potential_induced_C]
     type = Potential
     variable = potential
     concentration = conC
     dielectric =  138801 # 0.1388015 # inverse of dielectric constant, z_i  e  /epsilon
  [../]

  [./conA_diff]
    type = Diffusion
    variable = conA
  [../]

  [./conA_coupled]
    type = Concentration
    variable = conA
    potential = potential
    eqpotential = eqpot
    coefficient = -38.69565  # z_i  e  /KT
    bulkconc = 1
    gradchem = '0.0 0.0 0.0'
  [../]

  [./conC_diff]
    type = Diffusion
    variable = conC
  [../]

  [./conC_coupled]
    type = Concentration
    variable = conC
    potential = potential
    eqpotential = eqpot
    coefficient = 38.69565 # z_i  e  /KT
    bulkconc = 1
    gradchem = '0.0 0.0 0.0'
  [../]

[]

[AuxKernels]
  [./tot_conA]
    type = TotalConc
    variable = totConA
    eqpotential = eqpot
    concentration = conA
    coefficient = -38.69565
    bulkconc = 1
  [../]
  [./tot_conC]
    type = TotalConc
    variable = totConC
    eqpotential = eqpot
    concentration = conC
    coefficient =  38.69565
    bulkconc = 1
  [../]
[]

[Functions]
  [./bc_func]
     type = ParsedFunction
     value = 'sin(2*pi*t/100)*alpha+beta'
     vars = 'alpha beta'
     vals = '0.1 0.0'
  [../]

  [./charge_func]
     type = ParsedFunction
     value = '-sin(2*pi*t/100)*alpha'
     vars = 'alpha'
     vals = '230'
  [../]
[]


[BCs]
  [./eq_top]
    type = DirichletBC
    variable = eqpot
    boundary = top
    value    = 0.0
  [../]
  [./eq_bot]
    type = DirichletBC
    variable = eqpot
    boundary = bottom
    value    = 0.0
  [../]
  [./eq_gate_bc]
     type = NeumannBC
     #type = FunctionNeumannBC
     variable = eqpot
     boundary = 'right'
     value    = -230.0  # negative charged gate
     #function = charge_func
  [../]

  [./p_top]
    type = DirichletBC
    variable = potential
    boundary = top
    value    = -0.01
  [../]
  [./p_bottom]
    type = DirichletBC
    variable = potential
    boundary = bottom
    value    = 0.01
  [../]

  [./gate_bc]
     type = NeumannBC
     variable = potential
     boundary = 'right'
     value    = 0.0  # negative charged gate
  [../]

  [./inlet_conA]
    type = DirichletBC
    variable = conA
    boundary = bottom
    value = 0.0 # (C)
  [../]

  [./outlet_conA]
    type = DirichletBC
    variable = conA
    boundary = top
    value = 0.0 # (C)
  [../]

  [./bottomflow_conA]
    type = ConductionOutflow
    variable = conA
    potential = potential
    eqpotential = eqpot
    coefficient = -38.69565
    bulkconc = 1
    gradchem = '0.0 0.0 0.0'
    boundary = bottom
  [../]

  [./topflow_conA]
    type = ConductionOutflow
    variable = conA
    potential = potential
    eqpotential = eqpot
    coefficient = -38.69565
    bulkconc = 1
    gradchem = '0.0 0.0 0.0'
    boundary = top
  [../]

  [./inlet_conC]
    type = DirichletBC
    variable = conC
    boundary = bottom
    value = 0.0 # (C)
  [../]

  [./outlet_conC]
    type = DirichletBC
    variable = conC
    boundary = top
    value = 0.0 # (C)
  [../]

  [./bottomflow_conC]
    type = ConductionOutflow
    variable = conC
    potential = potential
    eqpotential = eqpot
    coefficient = 38.69565
    bulkconc = 1
    gradchem = '0.0 0.0 0.0'
    boundary = bottom
  [../]
  [./topflow_conC]
    type = ConductionOutflow
    variable = conC
    potential = potential
    eqpotential = eqpot
    coefficient = 38.69565
    bulkconc = 1
    gradchem = '0.0 0.0 0.0'
    boundary = top
  [../]

[]


#[UserObjects]
#  [./soln]
#    type = SolutionUserObject
#    system = nl0
#    mesh = pnp-restart_out_0005_mesh.xda
#    es = pnp-restart_out_0005.xda
#    system_variables = 'potential conA conC'
#    execute_on = initial
#  [../]
#[]

[Postprocessors]
  [./bottom_flux_ConA]
    type = SideFluxIntegralNp
    variable = conA
    potential = potential
    eqpotential = eqpot
    boundary = bottom
    coefficient = -38.69565
    bulkconc = 1
    gradchem = '0.0 0.0 0.0'
  [../]
  [./top_flux_ConA]
    type = SideFluxIntegralNp
    variable = conA
    potential = potential
    eqpotential = eqpot
    boundary = top
    coefficient = -38.69565
    bulkconc = 1
    gradchem = '0.0 0.0 0.0'
  [../]

  [./bottom_flux_ConC]
    type = SideFluxIntegralNp
    variable = conC
    potential = potential
    eqpotential = eqpot
    boundary = bottom
    coefficient = 38.69565
    bulkconc = 1
    gradchem = '0.0 0.0 0.0'
  [../]
  [./top_flux_ConC]
    type = SideFluxIntegralNp
    variable = conC
    potential = potential
    eqpotential = eqpot
    boundary = top
    coefficient = 38.69565
    bulkconc = 1
    gradchem = '0.0 0.0 0.0'
  [../]
[]

[Problem]
#  type = FEProblem
  coord_type = RZ
#  rz_coord_axis = Y
# restart_file_base = pnp-restart_out_cp/LATEST
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  #solve_type = PJFNK
  l_abs_step_tol = -1
  l_tol = 1.0e-5
  nl_rel_tol = 1.0e-3
  nl_max_its = 20
  l_max_its = 1000
  #petsc_options_iname = '-sub_pc_type'
  #petsc_options_value = '   hypre '
  #petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  #petsc_options_value = 'asm      31                  preonly       lu           8'

  #petsc_options_iname = '-ksp_gmres_restart '
  #petsc_options_value = '300                '
  #line_search = 'none'
[]

[Adaptivity]
  marker = 'errorfrac1 errorfrac2 errorfrac3 errorfrac4'
  steps = 2
  [./Indicators]
    [./error_eqp]
      type = GradientJumpIndicator
      variable = eqpot
    [../]
    [./error_p]
      type = GradientJumpIndicator
      variable = potential
    [../]
    [./error_conA]
      type = GradientJumpIndicator
      variable = conA
    [../]
    [./error_conC]
      type = GradientJumpIndicator
      variable = conC
    [../]
  [../]
  [./Markers]
    [./errorfrac1]
      type = ErrorFractionMarker
      refine = 0.5
      coarsen = 0
      indicator = error_eqp
    [../]
    [./errorfrac2]
      type = ErrorFractionMarker
      refine = 0.5
      coarsen = 0
      indicator = error_p
    [../]
    [./errorfrac3]
      type = ErrorFractionMarker
      refine = 0.5
      coarsen = 0
      indicator = error_conA
    [../]
    [./errorfrac4]
      type = ErrorFractionMarker
      refine = 0.5
      coarsen = 0
      indicator = error_conC
    [../]
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[Outputs]
  exodus = true
 # checkpoint = true
 # xda = true
[]
