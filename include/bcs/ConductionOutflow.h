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

#ifndef HEATCONDUCTIONOUTFLOW_H
#define HEATCONDUCTIONOUTFLOW_H

#include "IntegratedBC.h"


class ConductionOutflow;

template<>
InputParameters validParams<ConductionOutflow>();

/**
 * An IntegratedBC representing the "No BC" boundary condition for Heat Conduction.
 *
 * This is essentially -test*k*grad_u*normal... essentially completing the integration by parts.
 * This is a well accepted practice for truncating longer domains for convection/diffusion
 * problems as analyzed in: Griffiths, David F. "The ‘no boundary condition’outflow boundary condition." International journal for numerical methods in fluids 24.4 (1997): 393-411.
 */
class ConductionOutflow : public IntegratedBC
{
public:
  ConductionOutflow(const InputParameters & parameters);


protected:
  /// This is called to integrate the residual across the boundary
  virtual Real computeQpResidual();

  /// Optional (but recommended!) to compute the derivative of the residual with respect to _this_ variable
  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const VariableGradient & _potential_gradient;
  const VariableGradient & _eqpotential_gradient;
  const VariableValue & _eqpotential_value;
  unsigned int _potential_var;

  Real _coefficient;
  Real _bulkconc;
  RealVectorValue _gradchem;
};


#endif //HEATCONDUCTIONOUTFLOW_H
