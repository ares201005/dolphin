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

#include "Flux.h"

template<>
InputParameters validParams<Flux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Declare the options for a MooseEnum.
  // These options will be presented to the user in Peacock
  // and if something other than these options is in the input file
  // an error will be printed
  MooseEnum component("x y z");

  // Use the MooseEnum to add a parameter called "component"
  params.addRequiredParam<MooseEnum>("component", component, "The desired component of velocity.");

  // Add a "coupling paramater" to get a variable from the input file.
  params.addRequiredCoupledVar("potential",     "The poential field.");
  params.addRequiredCoupledVar("eqpotential", "The variable representing the pressure.");
  params.addRequiredCoupledVar("concentration", "The concentration profile");
  params.addRequiredParam<Real>("bulkconc", "The coefficient");
  params.addRequiredParam<Real>("coefficient", "The coefficient");
  params.addRequiredParam<RealVectorValue>("gradchem", "gradient of chemical potential");

  return params;
}

Flux::Flux(const InputParameters & parameters) :
    AuxKernel(parameters),

    // This will automatically convert the MooseEnum to an integer
    _component(getParam<MooseEnum>("component")),

    // Get the gradient of the variable
    _potential_gradient(coupledGradient("potential")),
    _eqpotential_gradient(coupledGradient("eqpotential")),
    _con_gradient(coupledGradient("concentration")),

    _con(coupledValue("concentration")),
    _eqpotential_value(coupledValue("eqpotential")),

    _bulkconc(getParam<Real>("bulkconc")),
    _coefficient(getParam<Real>("coefficient")),
    _gradchem(getParam<RealVectorValue>("gradchem"))

{
}

Real
Flux::computeValue()
{
  // Access the gradient of the pressure at this quadrature point
  // Then pull out the "component" of it we are looking for (x, y or z)
  // Note that getting a particular component of a gradient is done using the
  // parenthesis operator
  // Diff=
  return (_con_gradient[_qp](_component)+
          _coefficient*_con[_qp]*(_potential_gradient[_qp](_component)+_eqpotential_gradient[_qp](_component)) +
          _coefficient*_bulkconc*std::exp(-_coefficient*_eqpotential_value[_qp])*_potential_gradient[_qp](_component))*192.9448;
  // not finished
}
