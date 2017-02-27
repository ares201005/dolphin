/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "INSMomentumBaseForce.h"

template<>
InputParameters validParams<INSMomentumBaseForce>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  params.addRequiredCoupledVar("p", "pressure");

  params.addRequiredCoupledVar("eqpotential",     "The poential field.");
  params.addRequiredCoupledVar("potential", "the variable representing the potential");
  params.addRequiredCoupledVar("conA", "The variable representing the concentration of anion.");
  params.addRequiredCoupledVar("conC", "The variable representing the concentration of cation");

  // Required parameters
  params.addRequiredParam<Real>("coefficient", "The coefficient");
  params.addRequiredParam<Real>("bulkconc", "The bulk concentration");
  params.addRequiredParam<Real>("force", "force factor");

  params.addRequiredParam<Real>("mu", "dynamic viscosity");
  params.addRequiredParam<Real>("rho", "density");
  params.addRequiredParam<unsigned>("component", "0,1,2 depending on if we are solving the x,y,z component of the momentum equation");
  params.addParam<bool>("integrate_p_by_parts", true, "Allows simulations to be run with pressure BC if set to false");

  return params;
}



INSMomentumBaseForce::INSMomentumBaseForce(const InputParameters & parameters) :
  Kernel(parameters),

  // Coupled variables
  _u_vel(coupledValue("u")),
  _v_vel(coupledValue("v")),
  _w_vel(coupledValue("w")),
  _p(coupledValue("p")),

  _conA(coupledValue("conA")),
  _conC(coupledValue("conC")),
  _eq_pot(coupledValue("eqpotential")),

  // Gradients
  _grad_u_vel(coupledGradient("u")),
  _grad_v_vel(coupledGradient("v")),
  _grad_w_vel(coupledGradient("w")),
  _grad_p(coupledGradient("p")),

  _pot_grad(coupledGradient("potential")),
  _eq_pot_grad(coupledGradient("eqpotential")),

  // Variable numberings
  _u_vel_var_number(coupled("u")),
  _v_vel_var_number(coupled("v")),
  _w_vel_var_number(coupled("w")),
  _p_var_number(coupled("p")),

  _conA_var_number(coupled("conA")),
  _conC_var_number(coupled("conC")),
  _pot_var_number(coupled("potential")),
  _eq_pot_var_number(coupled("eqpotential")),

  // Required parameters
  _coefficient(getParam<Real>("coefficient")),
  _bulkconc(getParam<Real>("bulkconc")),
  _force(getParam<Real>("force")),

  _mu(getParam<Real>("mu")),
  _rho(getParam<Real>("rho")),
  _component(getParam<unsigned>("component")),
  _integrate_p_by_parts(getParam<bool>("integrate_p_by_parts"))

  // Material properties
  // _dynamic_viscosity(getMaterialProperty<Real>("dynamic_viscosity"))
{
}



Real INSMomentumBaseForce::computeQpResidual()
{
  // The convection part, rho * (u.grad) * u_component * v.
  // Note: _grad_u is the gradient of the _component entry of the velocity vector.
  Real convective_part = _rho *
    (_u_vel[_qp]*_grad_u[_qp](0) +
     _v_vel[_qp]*_grad_u[_qp](1) +
     _w_vel[_qp]*_grad_u[_qp](2)) * _test[_i][_qp];

  // The pressure part, -p (div v) or (dp/dx_{component}) * test if not integrated by parts.
  Real pressure_part = 0.;
  if (_integrate_p_by_parts)
    pressure_part = -_p[_qp] * _grad_test[_i][_qp](_component);
  else
    pressure_part = _grad_p[_qp](_component) * _test[_i][_qp];

  // Call derived class to compute the viscous contribution.
  Real viscous_part = computeQpResidualViscousPart();

  // Body force term.  For truly incompressible flow, this term is constant, and
  // since it is proportional to g, can be written as the gradient of some scalar
  // and absorbed into the pressure definition.
  // Real body_force_part = - _rho * _gravity(_component);

  Real body_force_part = _force * ( _pot_grad[_qp](_component) + _eq_pot_grad[_qp](_component) ) * 
                         (_conC[_qp] + _bulkconc*std::exp(-_coefficient*_eq_pot[_qp]) - 
                          _conA[_qp] - _bulkconc*std::exp( _coefficient*_eq_pot[_qp]))*_test[_i][_qp];

  convective_part = 0;

  return convective_part + pressure_part + viscous_part + body_force_part;
}




Real INSMomentumBaseForce::computeQpJacobian()
{
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);

  // Convective part
  Real convective_part = _rho * ((U*_grad_phi[_j][_qp]) + _phi[_j][_qp]*_grad_u[_qp](_component)) * _test[_i][_qp];

  // Call derived class to compute the viscous contribution.
  Real viscous_part = computeQpJacobianViscousPart();

  convective_part = 0;

  return convective_part + viscous_part;
}




Real INSMomentumBaseForce::computeQpOffDiagJacobian(unsigned jvar)
{
  // In Stokes/Laplacian version, off-diag Jacobian entries wrt u,v,w are zero
  if (jvar == _u_vel_var_number)
  {
    Real convective_part = _phi[_j][_qp] * _grad_u[_qp](0) * _test[_i][_qp];

    // Call derived class to compute the viscous contribution.
    Real viscous_part = computeQpOffDiagJacobianViscousPart(jvar);

    convective_part = 0;

    return convective_part + viscous_part;
  }

  else if (jvar == _v_vel_var_number)
  {
    Real convective_part = _phi[_j][_qp] * _grad_u[_qp](1) * _test[_i][_qp];

    // Call derived class to compute the viscous contribution.
    Real viscous_part = computeQpOffDiagJacobianViscousPart(jvar);

    convective_part = 0;

    return convective_part + viscous_part;
  }

  else if (jvar == _w_vel_var_number)
  {
    Real convective_part = _phi[_j][_qp] * _grad_u[_qp](2) * _test[_i][_qp];

    // Call derived class to compute the viscous contribution.
    Real viscous_part = computeQpOffDiagJacobianViscousPart(jvar);

    convective_part = 0;

    return convective_part + viscous_part;
  }

  else if (jvar == _p_var_number)
  {
    if (_integrate_p_by_parts)
      return -_phi[_j][_qp] * _grad_test[_i][_qp](_component);
    else
      return _grad_phi[_j][_qp](_component) * _test[_i][_qp];
  }

  else if (jvar == _conA_var_number)
    return -_force * _phi[_j][_qp] * (_pot_grad[_qp](_component) + _eq_pot_grad[_qp](_component)) * _test[_i][_qp];

  else if (jvar == _conC_var_number)
    return _force * _phi[_j][_qp] * (_pot_grad[_qp](_component) + _eq_pot_grad[_qp](_component)) * _test[_i][_qp];

  else if (jvar == _pot_var_number)
    return _force * (_conC[_qp] + _bulkconc*std::exp(-_coefficient*_eq_pot[_qp]) -
                     _conA[_qp] - _bulkconc*std::exp( _coefficient*_eq_pot[_qp])) * _grad_phi[_j][_qp](_component) * _test[_i][_qp];

  else if (jvar == _eq_pot_var_number)
  {
    Real partA =  _force * (_conC[_qp] + _bulkconc*std::exp(-_coefficient*_eq_pot[_qp]) -
                            _conA[_qp] - _bulkconc*std::exp( _coefficient*_eq_pot[_qp])) * _grad_phi[_j][_qp](_component) * _test[_i][_qp];
 
    Real partB = -_force * ( _pot_grad[_qp](_component) + _eq_pot_grad[_qp](_component) ) * 
                 _coefficient * _bulkconc * (std::exp( _coefficient*_eq_pot[_qp]) + std::exp(-_coefficient*_eq_pot[_qp])) *
                 _phi[_j][_qp] * _test[_i][_qp];

    return partA + partB;
  }

  else
    return 0;
}
