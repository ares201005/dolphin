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

#ifndef SIDEFLUXINTEGRALNP_H
#define SIDEFLUXINTEGRALNP_H

// MOOSE includes
#include "SideIntegralVariablePostprocessor.h"

// Forward Declarations
class SideFluxIntegralNp;

template<>
InputParameters validParams<SideFluxIntegralNp>();

/**
 * This postprocessor computes a side integral of the mass flux.
 */
class SideFluxIntegralNp : public SideIntegralVariablePostprocessor
{
public:
  SideFluxIntegralNp(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

  const VariableGradient & _potential_gradient;
  const VariableGradient & _eqpotential_gradient;
  const VariableValue & _eqpotential_value;

  Real _coefficient;
  Real _bulkconc;
  RealVectorValue _gradchem;

};

#endif // SIDEFLUXINTEGRAL_H
