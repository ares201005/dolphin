/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMFORCELAPLACEFORMRZ_H
#define INSMOMENTUMFORCELAPLACEFORMRZ_H

#include "INSMomentumForceLaplaceForm.h"

// Forward Declarations
class INSMomentumForceLaplaceFormRZ;

template<>
InputParameters validParams<INSMomentumForceLaplaceFormRZ>();

/**
 * This class computes additional momentum equation residual and
 * Jacobian contributions for the incompressible Navier-Stokes
 * momentum equation in RZ (axisymmetric cylindrical) coordinates,
 * using the "Laplace" form of the governing equations.
 */
class INSMomentumForceLaplaceFormRZ : public INSMomentumForceLaplaceForm
{
public:
  INSMomentumForceLaplaceFormRZ(const InputParameters & parameters);

  virtual ~INSMomentumForceLaplaceFormRZ(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);
};


#endif
