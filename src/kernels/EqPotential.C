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

#include "EqPotential.h"


template<>
InputParameters validParams<EqPotential>()
{
  // Start with the parameters from our parent
  InputParameters params = validParams<Kernel>();

  params.addRequiredParam<Real>("dielectric", "The dielectric constant");
  params.addRequiredParam<Real>("coefficient", "The coefficient");
  params.addRequiredParam<Real>("bulkconc", "The bulk concentration");

  return params;
}


EqPotential::EqPotential(const InputParameters & parameters) :
    Kernel(parameters),

    _dielectric(getParam<Real>("dielectric")),
    _coefficient(getParam<Real>("coefficient")),
    _bulkconc(getParam<Real>("bulkconc"))
{
}

Real
EqPotential::computeQpResidual()
{
  // Use the MaterialProperty references we stored earlier
  return _dielectric *_bulkconc* (std::exp(_coefficient*_u[_qp])-std::exp(-_coefficient*_u[_qp])) *_test[_i][_qp];
}

Real
EqPotential::computeQpJacobian()
{
  return _dielectric*_bulkconc*(std::exp(_coefficient*_u[_qp])+std::exp(-_coefficient*_u[_qp])) * _coefficient*_phi[_j][_qp] *_test[_i][_qp];
}
