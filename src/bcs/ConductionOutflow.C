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

#include "ConductionOutflow.h"

template<>
InputParameters validParams<ConductionOutflow>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addRequiredCoupledVar("potential", "The variable representing the pressure.");
  params.addRequiredCoupledVar("eqpotential", "The variable representing the pressure.");
  params.addRequiredParam<Real>("coefficient", "The coefficient");
  params.addRequiredParam<Real>("bulkconc", "The bulk concentration");
  params.addRequiredParam<RealVectorValue>("gradchem", "gradient of chemical potential");

  return params;
}

ConductionOutflow::ConductionOutflow(const InputParameters & parameters) :
  IntegratedBC(parameters),

  _potential_gradient(coupledGradient("potential")),
  _eqpotential_gradient(coupledGradient("eqpotential")),
  _eqpotential_value(coupledValue("eqpotential")),
  _coefficient(getParam<Real>("coefficient")),
  _bulkconc(getParam<Real>("bulkconc")),
  _gradchem(getParam<RealVectorValue>("gradchem"))
  // IntegratedBCs can retrieve material properties!
{}

Real
ConductionOutflow::computeQpResidual()
{
  // See: Griffiths, David F. "The ‘no boundary condition’outflow boundary condition." International journal for numerical methods in fluids 24.4 (1997): 393-411.
  return -_test[_i][_qp]*(_grad_u[_qp]+_coefficient*_u[_qp]*(_potential_gradient[_qp]+_eqpotential_gradient[_qp]+_gradchem)+
         +_coefficient*_bulkconc*std::exp(-_coefficient*_eqpotential_value[_qp])*(_potential_gradient[_qp]+_gradchem))*_normals[_qp];
}

Real
ConductionOutflow::computeQpJacobian()
{
  // Derivative of the residual with respect to "u"
  return -_test[_i][_qp]*(_grad_phi[_j][_qp]+_coefficient*_phi[_j][_qp]*(_potential_gradient[_qp]+_eqpotential_gradient[_qp]+_gradchem))*_normals[_qp];
}

Real
ConductionOutflow::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Derivative of the residual with respect to "u"
  if (jvar == _potential_var)
  {
    return -_test[_i][_qp]*_coefficient*(_u[_qp]+_bulkconc*std::exp(-_coefficient*_eqpotential_value[_qp]))*_grad_phi[_j][_qp]*_normals[_qp];
  }
  return 0.0;
}
