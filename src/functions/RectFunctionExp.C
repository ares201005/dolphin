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

#include "RectFunctionExp.h"

template<>
InputParameters validParams<RectFunctionExp>()
{
  InputParameters params = validParams<Function>();
  params.addParam<Real>("alpha", 1.0, "The value of alpha");
  params.addParam<Real>("tau", 1.0, "The value of tau");
  return params;
}

RectFunctionExp::RectFunctionExp(const InputParameters & parameters) :
    Function(parameters),
    _alpha(getParam<Real>("alpha")),
    _tau(getParam<Real>("tau"))
{}

Real
RectFunctionExp::value(Real t, const Point & /*p*/)
{
  Real t0 = int(t/_tau) * _tau;
  Real t1 = (int(t/_tau) +1) * _tau;
  Real at = _tau/25.0;

  if( int(t/_tau) % 2 == 0 ) {
    return _alpha *(1.0-std::exp(-(t-t0)/at))*(1.0-exp(-(t1-t)/at));
  }
  else {
    return -_alpha *(1.0-std::exp(-(t-t0)/at))*(1.0-exp(-(t1-t)/at));
  }
}
