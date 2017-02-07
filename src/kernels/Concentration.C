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

#include "Concentration.h"

template<>
InputParameters validParams<Concentration>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredCoupledVar("potential", "The variable representing the pressure.");
  params.addRequiredCoupledVar("eqpotential", "The variable representing the pressure.");
  params.addRequiredParam<Real>("coefficient", "The coefficient");
  params.addRequiredParam<Real>("bulkconc", "The bulk concentrantion");
  params.addRequiredParam<RealVectorValue>("gradchem", "gradient of chemical potential");

  return params;
}

Concentration::Concentration(const InputParameters & parameters) :
    Kernel(parameters),

    // Couple to the gradient of the pressure
    _potential_gradient(coupledGradient("potential")),
    _eqpotential_gradient(coupledGradient("eqpotential")),
    _eqpotential_value(coupledValue("eqpotential")),

    // Save off the coupled variable identifier for use in
    // computeQpOffDiagJacobian
    _potential_var(coupled("potential")),
    _eqpotential_var(coupled("eqpotential")),

    _coefficient(getParam<Real>("coefficient")),
    _bulkconc(getParam<Real>("bulkconc")),

    _gradchem(getParam<RealVectorValue>("gradchem"))

    // Grab necessary material properties
{
}

Real
Concentration::computeQpResidual()
{
  // See also: E. Majchrzak and L. Turchan, "The Finite Difference
  // Method For Transient Convection Diffusion", Scientific Research
  // of the Institute of Mathematics and Computer Science, vol. 1,
  // no. 11, 2012, pp. 63-72.
  // http://srimcs.im.pcz.pl/2012_1/art_07.pdf

  return _coefficient * (_u[_qp] * (_potential_gradient[_qp]+_eqpotential_gradient[_qp]+_gradchem)+
         _bulkconc*std::exp(-_coefficient*_eqpotential_value[_qp]) * (_potential_gradient[_qp]+_gradchem)) * _grad_test[_i][_qp];
}

Real
Concentration::computeQpJacobian()
{
  return _coefficient * _phi[_j][_qp] * (_potential_gradient[_qp] + _eqpotential_gradient[_qp]+_gradchem)* _grad_test[_i][_qp];
}

Real
Concentration::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _potential_var)
  {
    return _coefficient * (_u[_qp] + _bulkconc*std::exp(-_coefficient*_eqpotential_value[_qp])) * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
  }
  //else if (jvar == _eqpotential_var)
  //{
  //  RealVectorValue tmp_grad1 =  -std::exp(-_coefficient*_eqpotential_value[_qp])*_coefficient*_phi[_j][_qp] *_potential_gradient[_qp]; 
  //  RealVectorValue tmp_grad2 = _u[_qp] * _grad_phi[_j][_qp];
  //  return _coefficient * (tmp_grad1 + tmp_grad2) * _grad_test[_i][_qp];
  //}
  return 0.0;
}
