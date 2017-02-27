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

#ifndef HEATCONDUCTIONOUTFLOWNS_H
#define HEATCONDUCTIONOUTFLOWNS_H

#include "IntegratedBC.h"


class ConductionOutflowNs;

template<>
InputParameters validParams<ConductionOutflowNs>();

/**
 * An IntegratedBC representing the "No BC" boundary condition for Heat Conduction.
 *
 * This is essentially -test*k*grad_u*normal... essentially completing the integration by parts.
 * This is a well accepted practice for truncating longer domains for convection/diffusion
 * problems as analyzed in: Griffiths, David F. "The ‘no boundary condition’outflow boundary condition." International journal for numerical methods in fluids 24.4 (1997): 393-411.
 */
class ConductionOutflowNs : public IntegratedBC
{
public:
  ConductionOutflowNs(const InputParameters & parameters);


protected:
  /// This is called to integrate the residual across the boundary
  virtual Real computeQpResidual();

  /// Optional (but recommended!) to compute the derivative of the residual with respect to _this_ variable
  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const VariableGradient & _potential_gradient;
  const VariableGradient & _eqpotential_gradient;
  const VariableValue & _eqpotential_value;

  Real _bulkconc;
  Real _coefficient;
  RealVectorValue _gradchem;

  Real _over_diff;

  // Coupled variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;

  // Variable numberings
  unsigned int _potential_var;

  unsigned _u_vel_var_number;
  unsigned _v_vel_var_number;
  unsigned _w_vel_var_number;
};


#endif //HEATCONDUCTIONOUTFLOWNS_H
