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

#ifndef SIDEFLUXINTEGRALNS_H
#define SIDEFLUXINTEGRALNS_H

// MOOSE includes
#include "SideIntegralVariablePostprocessor.h"

// Forward Declarations
class SideFluxIntegralNs;

template<>
InputParameters validParams<SideFluxIntegralNs>();

/**
 * This postprocessor computes a side integral of the mass flux.
 */
class SideFluxIntegralNs : public SideIntegralVariablePostprocessor
{
public:
  SideFluxIntegralNs(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

  const VariableGradient & _potential_gradient;
  const VariableGradient & _eqpotential_gradient;
  const VariableValue & _eqpotential_value;

  Real _coefficient;
  Real _bulkconc;
  RealVectorValue _gradchem;

  Real _over_diff;

  // Coupled variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
};

#endif // SIDEFLUXINTEGRAL_H
