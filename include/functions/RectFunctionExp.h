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

#ifndef RECTFUNCTIONEXP_H
#define RECTFUNCTIONEXP_H

#include "Function.h"

class RectFunctionExp;

template<>
InputParameters validParams<RectFunctionExp>();

class RectFunctionExp : public Function
{
public:
  RectFunctionExp(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p);

protected:
  Real _alpha;
  Real _tau;
};

#endif //RECTFUNCTIONEXP_H
