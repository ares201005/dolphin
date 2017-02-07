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

#include "SideFluxIntegralNp.h"

template<>
InputParameters validParams<SideFluxIntegralNp>()
{
  InputParameters params = validParams<SideIntegralVariablePostprocessor>();
  params.addRequiredCoupledVar("potential", "The variable representing the pressure.");
  params.addRequiredCoupledVar("eqpotential", "The variable representing the pressure.");
  params.addRequiredParam<Real>("coefficient", "The coefficient");
  params.addRequiredParam<Real>("bulkconc", "The coefficient");
  params.addRequiredParam<RealVectorValue>("gradchem", "gradient of chemical potential");
  return params;
}

SideFluxIntegralNp::SideFluxIntegralNp(const InputParameters & parameters) :
    SideIntegralVariablePostprocessor(parameters),
  _potential_gradient(coupledGradient("potential")),
  _eqpotential_gradient(coupledGradient("eqpotential")),
  _eqpotential_value(coupledValue("eqpotential")),
  _coefficient(getParam<Real>("coefficient")),
  _bulkconc(getParam<Real>("bulkconc")),
  _gradchem(getParam<RealVectorValue>("gradchem"))
{}

Real
SideFluxIntegralNp::computeQpIntegral()
{

  return -(_grad_u[_qp]+_coefficient*_u[_qp]*(_potential_gradient[_qp]+_eqpotential_gradient[_qp]+_gradchem)+
  +_coefficient*_bulkconc*std::exp(-_coefficient*_eqpotential_value[_qp])*(_potential_gradient[_qp]+_gradchem))*_normals[_qp];
}

