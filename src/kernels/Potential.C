/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "Potential.h"


template<>
InputParameters validParams<Potential>()
{
  // Start with the parameters from our parent
  InputParameters params = validParams<Kernel>();

  params.addRequiredCoupledVar("concentration", "The variable representing the concentration of anion or canion.");
  // No parameters are necessary here because we're going to get
  // permeability and viscosity from the Material
  // so we just return params...
  params.addRequiredParam<Real>("dielectric", "The dielectric constant");
  params.addRequiredParam<Real>("charge", "The charge");

  return params;
}


Potential::Potential(const InputParameters & parameters) :
    Kernel(parameters),

    _concentration_value(coupledValue("concentration")),
    _concentration_var(coupled("concentration")),
    // Get the permeability and viscosity from the Material system
    // This returns a MaterialProperty<Real> reference that we store
    // in the class and then index into in computeQpResidual/Jacobian....
  //  _permeability(getMaterialProperty<Real>("permeability"))
    _dielectric(getParam<Real>("dielectric")),
    _charge(getParam<Real>("charge"))
{
}

Real
Potential::computeQpResidual()
{
  // Use the MaterialProperty references we stored earlier
 // return -_permeability[_qp] * _concentration_value[_qp]*_test[_i][_qp];
  return -_charge*_dielectric * _concentration_value[_qp]*_test[_i][_qp];
}

Real
Potential::computeQpJacobian()
{
  return 0.0;  // Use the MaterialProperty references we stored earlier
}

Real
Potential::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Use the MaterialProperty references we stored earlier
  if (jvar == _concentration_var)
  {
   // return -_permeability[_qp] * _phi[_j][_qp]*_test[_i][_qp];
    return -_charge*_dielectric * _phi[_j][_qp]*_test[_i][_qp];
  }
  return 0.0;
}
