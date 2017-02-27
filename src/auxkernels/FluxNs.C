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

#include "FluxNs.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<FluxNs>()
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
  params.addRequiredParam<Real>("overdiff", "inverse of diffusion coefficient");

  // Coupled variables
  //params.addCoupledVar("u", "x-velocity");
  params.addCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", "z-velocity"); // only required in 3D


  return params;
}

FluxNs::FluxNs(const InputParameters & parameters) :
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
    _gradchem(getParam<RealVectorValue>("gradchem")),
    _over_diff(getParam<Real>("overdiff")),

  // Coupled variables
    _u_vel(coupledValue("u")),
    _v_vel(_mesh.dimension() >= 2 ? coupledValue("v") : _zero),
    _w_vel(_mesh.dimension() == 3 ? coupledValue("w") : _zero)

{
}

Real
FluxNs::computeValue()
{
  // Access the gradient of the pressure at this quadrature point
  // Then pull out the "component" of it we are looking for (x, y or z)
  // Note that getting a particular component of a gradient is done using the
  // parenthesis operator
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  //RealVectorValue U(0, 0, 0);

  return _con_gradient[_qp](_component)+_coefficient*_con[_qp]*_potential_gradient[_qp](_component) - _over_diff*_con[_qp]*U(_component);
  // not finished
}
