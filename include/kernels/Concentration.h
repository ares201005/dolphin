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

#ifndef CONCENTRATION_H
#define CONCENTRATION_H

#include "Kernel.h"

// Forward Declaration
class Concentration;

template<>
InputParameters validParams<Concentration>();

/**
 * Kernel which implements the convective term in the transient heat
 * conduction equation, and provides coupling with the Darcy pressure
 * equation.
 */
class Concentration : public Kernel
{
public:
  Concentration(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// The gradient of pressure
  const VariableGradient & _potential_gradient;
  const VariableGradient & _eqpotential_gradient;
  const VariableValue & _eqpotential_value;

  /// Coupling identifier for the pressure.  This is used to uniquely
  /// identify a coupled variable
  unsigned int _potential_var;
  unsigned int _eqpotential_var;

  Real _coefficient;
  Real _bulkconc;

  RealVectorValue _gradchem;

  /**
   * These references will be set by the initialization list so that
   * values can be pulled from the Material system.
   */
};

#endif //CONCENTRATIONA_H
