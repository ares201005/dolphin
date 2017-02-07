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

#ifndef POTENTIAL_H
#define POTENTIAL_H

// Including the "Diffusion" Kernel here so we can extend it
#include "Kernel.h"

class Potential;

template<>
InputParameters validParams<Potential>();


/**
 * Represents K/mu * grad_u * grad_phi
 *
 * We are inheriting from Diffusion instead of from Kernel because
 * the grad_u * grad_phi is already coded in there and all we
 * need to do is specialize that calculation by multiplying by K/mu
 */
class Potential : public Kernel
{
public:
  Potential(const InputParameters & parameters);

protected:
  /**
   * Kernels _must_ override computeQpResidual()
   */
  virtual Real computeQpResidual();

  /**
   * This is optional (but recommended!)
   */
  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const VariableValue & _concentration_value;

  unsigned int _concentration_var;

  Real _dielectric;
  /**
   * These references will be set by the initialization list so that
   * values can be pulled from the Material system.
   */
  //const MaterialProperty<Real> & _permeability;
};


#endif /* POTENTIAL_H */