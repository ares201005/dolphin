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

#ifndef TOTALCONC_H
#define TOTALCONC_H

#include "AuxKernel.h"

//Forward Declarations
class TotalConc;

template<>
InputParameters validParams<TotalConc>();

/**
 * Constant auxiliary value
 */
class TotalConc : public AuxKernel
{
public:
  TotalConc(const InputParameters & parameters);

  virtual ~TotalConc() {}

protected:
  /**
   * AuxKernels MUST override computeValue.  computeValue() is called on
   * every quadrature point.  For Nodal Auxiliary variables those quadrature
   * points coincide with the nodes.
   */
  virtual Real computeValue();

  /// Will hold 0,1,2 for x,y,z
  /// The gradient of a coupled variable
  const VariableValue & _eqpotential;
  const VariableValue & _con;
  Real _coefficient;
  Real _bulkconc;

  
  /// Holds the permeability and viscosity from the material system
};

#endif  //TOTALCONC_H
