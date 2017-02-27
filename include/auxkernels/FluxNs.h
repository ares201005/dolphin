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

#ifndef FLUXNS_H
#define FLUXNS_H

#include "AuxKernel.h"

//Forward Declarations
class FluxNs;

template<>
InputParameters validParams<FluxNs>();

/**
 * Constant auxiliary value
 */
class FluxNs : public AuxKernel
{
public:
  FluxNs(const InputParameters & parameters);

  virtual ~FluxNs() {}

protected:
  /**
   * AuxKernels MUST override computeValue.  computeValue() is called on
   * every quadrature point.  For Nodal Auxiliary variables those quadrature
   * points coincide with the nodes.
   */
  virtual Real computeValue();

  /// Will hold 0,1,2 for x,y,z
  int _component;

  /// The gradient of a coupled variable
  const VariableGradient & _potential_gradient;
  const VariableGradient & _eqpotential_gradient;
  const VariableGradient & _con_gradient;

  const VariableValue & _con;
  const VariableValue & _eqpotential_value;

  Real _bulkconc;
  Real _coefficient;
  RealVectorValue _gradchem;
  Real _over_diff;

  // Coupled variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;

  
  /// Holds the permeability and viscosity from the material system
};

#endif //DARCYVELOCITY_H
