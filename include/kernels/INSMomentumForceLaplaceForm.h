/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMFORCELAPLACEFORM_H
#define INSMOMENTUMFORCELAPLACEFORM_H

#include "INSMomentumBaseForce.h"

// Forward Declarations
class INSMomentumForceLaplaceForm;

template<>
InputParameters validParams<INSMomentumForceLaplaceForm>();

/**
 * This class computes momentum equation residual and Jacobian viscous
 * contributions for the "Laplacian" form of the governing equations.
 */
class INSMomentumForceLaplaceForm : public INSMomentumBaseForce
{
public:
  INSMomentumForceLaplaceForm(const InputParameters & parameters);

  virtual ~INSMomentumForceLaplaceForm(){}

protected:
  virtual Real computeQpResidualViscousPart() override;
  virtual Real computeQpJacobianViscousPart() override;
  virtual Real computeQpOffDiagJacobianViscousPart(unsigned jvar) override;
};


#endif
