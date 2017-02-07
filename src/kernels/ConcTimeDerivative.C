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

#include "ConcTimeDerivative.h"

template<>
InputParameters validParams<ConcTimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  params.addParam<Real>("time_coefficient", 1.0, "Time Coefficient");
  return params;
}

ConcTimeDerivative::ConcTimeDerivative(const InputParameters & parameters) :
    TimeDerivative(parameters),
    // This kernel expects an input parameter named "time_coefficient"
    _time_coefficient(getParam<Real>("time_coefficient"))
{}

Real
ConcTimeDerivative::computeQpResidual()
{
  return _time_coefficient*TimeDerivative::computeQpResidual();
}

Real
ConcTimeDerivative::computeQpJacobian()
{
  return _time_coefficient*TimeDerivative::computeQpJacobian();
}
