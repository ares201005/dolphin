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

#ifndef RECTFUNCTION_H
#define RECTFUNCTION_H

#include "Function.h"

class RectFunction;

template<>
InputParameters validParams<RectFunction>();

class RectFunction : public Function
{
public:
  RectFunction(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p);

protected:
  Real _alpha;
  Real _tau;
};

#endif //RECTFUNCTION_H
