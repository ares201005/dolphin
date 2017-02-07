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

#include "RectFunction.h"

template<>
InputParameters validParams<RectFunction>()
{
  InputParameters params = validParams<Function>();
  params.addParam<Real>("alpha", 1.0, "The value of alpha");
  params.addParam<Real>("tau", 1.0, "The value of tau");
  return params;
}

RectFunction::RectFunction(const InputParameters & parameters) :
    Function(parameters),
    _alpha(getParam<Real>("alpha")),
    _tau(getParam<Real>("tau"))
{}

Real
RectFunction::value(Real t, const Point & /*p*/)
{
  if( int(t/_tau) % 2 == 0 ) {
    return _alpha;
  }
  else {
    return _alpha*(-1.0);
  }
}
