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

#include "ConvectionVel.h"
#include "MooseMesh.h"

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<ConvectionVel>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredParam<Real>("overdiff", "inverse of diffusion coefficient");
  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  //params.addCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", "z-velocity"); // only required in 3D


  return params;
}

ConvectionVel::ConvectionVel(const InputParameters & parameters) :
  Kernel(parameters),

  _over_diff(getParam<Real>("overdiff")),
  // Coupled variables
  _u_vel(coupledValue("u")),
  _v_vel(_mesh.dimension() >= 2 ? coupledValue("v") : _zero),
  _w_vel(_mesh.dimension() == 3 ? coupledValue("w") : _zero),

  // Variable numberings
  _u_vel_var_number(coupled("u")),
  _v_vel_var_number(_mesh.dimension() >= 2 ? coupled("v") : libMesh::invalid_uint),
  _w_vel_var_number(_mesh.dimension() == 3 ? coupled("w") : libMesh::invalid_uint)

{
}

Real ConvectionVel::computeQpResidual()
{
  // velocity * _grad_u[_qp] is actually doing a dot product

  RealVectorValue U1(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);

  return -_over_diff * _grad_test[_i][_qp] * (U1 * _u[_qp]);
}

Real ConvectionVel::computeQpJacobian()
{
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);

  // the partial derivative of _grad_u is just _grad_phi[_j]
  return -_over_diff * _grad_test[_i][_qp]*(U * _phi[_j][_qp]);
}


Real ConvectionVel::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _u_vel_var_number)
  {
    return -_over_diff * _grad_test[_i][_qp](0) * _phi[_j][_qp] * _u[_qp];
  }
  
  else if (jvar == _v_vel_var_number)
  {
    return -_over_diff * _grad_test[_i][_qp](1) * _phi[_j][_qp] * _u[_qp];
  }
  
  else if (jvar == _w_vel_var_number)
  {
    return -_over_diff * _grad_test[_i][_qp](2) * _phi[_j][_qp] * _u[_qp];
  }
  else
    return 0;
}
