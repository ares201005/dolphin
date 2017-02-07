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

#include "TotalConc.h"

template<>
InputParameters validParams<TotalConc>()
{
  InputParameters params = validParams<AuxKernel>();

  // Declare the options for a MooseEnum.
  // These options will be presented to the user in Peacock
  // and if something other than these options is in the input file
  // an error will be printed

  // Add a "coupling paramater" to get a variable from the input file.
  params.addRequiredCoupledVar("eqpotential",     "The poential field.");
  params.addRequiredCoupledVar("concentration", "The concentration profile");
  params.addRequiredParam<Real>("coefficient", "The coefficient");
  params.addRequiredParam<Real>("bulkconc", "The coefficient");

  return params;
}

TotalConc::TotalConc(const InputParameters & parameters) :
    AuxKernel(parameters),

    // This will automatically convert the MooseEnum to an integer

    // Get the gradient of the variable
    _eqpotential(coupledValue("eqpotential")),
    _con(coupledValue("concentration")),
    _coefficient(getParam<Real>("coefficient")),
    _bulkconc(getParam<Real>("bulkconc"))

{
}

Real
TotalConc::computeValue()
{
  // Access the gradient of the pressure at this quadrature point
  // Then pull out the "component" of it we are looking for (x, y or z)
  // Note that getting a particular component of a gradient is done using the
  // parenthesis operator
  return _con[_qp]+_bulkconc*std::exp(-_coefficient*_eqpotential[_qp]);
}
