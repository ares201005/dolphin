/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "INSMomentumForceLaplaceForm.h"

template<>
InputParameters validParams<INSMomentumForceLaplaceForm>()
{
  InputParameters params = validParams<INSMomentumBaseForce>();
  return params;
}



INSMomentumForceLaplaceForm::INSMomentumForceLaplaceForm(const InputParameters & parameters) :
  INSMomentumBaseForce(parameters)
{
}



Real INSMomentumForceLaplaceForm::computeQpResidualViscousPart()
{
  // Simplified version: mu * Laplacian(u_component)
  return _mu * (_grad_u[_qp] * _grad_test[_i][_qp]);
}



Real INSMomentumForceLaplaceForm::computeQpJacobianViscousPart()
{
  // Viscous part, Laplacian version
  return _mu * (_grad_phi[_j][_qp] * _grad_test[_i][_qp]);
}



Real INSMomentumForceLaplaceForm::computeQpOffDiagJacobianViscousPart(unsigned /*jvar*/)
{
  return 0.;
}
