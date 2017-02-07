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

#ifndef DARCYVELOCITY_H
#define DARCYVELOCITY_H

#include "AuxKernel.h"

//Forward Declarations
class Flux;

template<>
InputParameters validParams<Flux>();

/**
 * Constant auxiliary value
 */
class Flux : public AuxKernel
{
public:
  Flux(const InputParameters & parameters);

  virtual ~Flux() {}

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
  const VariableGradient & _con_gradient;
  const VariableValue & _con;
  Real _coefficient;
  RealVectorValue _gradchem;

  
  /// Holds the permeability and viscosity from the material system
};

#endif //DARCYVELOCITY_H
