/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMFORCETRACTIONFORM_H
#define INSMOMENTUMFORCETRACTIONFORM_H

#include "INSMomentumBaseForce.h"

// Forward Declarations
class INSMomentumForceTractionForm;

template<>
InputParameters validParams<INSMomentumForceTractionForm>();

/**
 * This class computes momentum equation residual and Jacobian viscous
 * contributions for the "traction" form of the governing equations.
 */
class INSMomentumForceTractionForm : public INSMomentumBaseForce
{
public:
  INSMomentumForceTractionForm(const InputParameters & parameters);

  virtual ~INSMomentumForceTractionForm(){}

protected:
  virtual Real computeQpResidualViscousPart() override;
  virtual Real computeQpJacobianViscousPart() override;
  virtual Real computeQpOffDiagJacobianViscousPart(unsigned jvar) override;
};


#endif
