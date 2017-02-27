/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMFORCETRACTIONFORMRZ_H
#define INSMOMENTUMFORCETRACTIONFORMRZ_H

#include "INSMomentumForceTractionForm.h"

// Forward Declarations
class INSMomentumForceTractionFormRZ;

template<>
InputParameters validParams<INSMomentumForceTractionFormRZ>();

/**
 * This class computes additional momentum equation residual and
 * Jacobian contributions for the incompressible Navier-Stokes
 * momentum equation in RZ (axisymmetric cylindrical) coordinates.
 */
class INSMomentumForceTractionFormRZ : public INSMomentumForceTractionForm
{
public:
  INSMomentumForceTractionFormRZ(const InputParameters & parameters);

  virtual ~INSMomentumForceTractionFormRZ(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);
};


#endif
