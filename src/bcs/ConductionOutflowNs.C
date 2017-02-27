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

#include "ConductionOutflowNs.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<ConductionOutflowNs>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addRequiredCoupledVar("potential", "The variable representing the pressure.");
  params.addRequiredCoupledVar("eqpotential", "The variable representing the pressure.");
  params.addRequiredParam<Real>("coefficient", "The coefficient");
  params.addRequiredParam<Real>("bulkconc", "The bulk concentration");
  params.addRequiredParam<RealVectorValue>("gradchem", "gradient of chemical potential");

  params.addRequiredParam<Real>("overdiff", "inverse of diffusion coefficient");

  // Coupled variables
  //params.addCoupledVar("u", "x-velocity");
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", "z-velocity"); // only required in 3D


  return params;
}

ConductionOutflowNs::ConductionOutflowNs(const InputParameters & parameters) :
  IntegratedBC(parameters),

  _potential_gradient(coupledGradient("potential")),
  _eqpotential_gradient(coupledGradient("eqpotential")),
  _eqpotential_value(coupledValue("eqpotential")),
  _bulkconc(getParam<Real>("bulkconc")),
  _coefficient(getParam<Real>("coefficient")),
  _gradchem(getParam<RealVectorValue>("gradchem")),
  _over_diff(getParam<Real>("overdiff")),
  // IntegratedBCs can retrieve material properties!

  // Coupled variables
   _u_vel(coupledValue("u")),
   _v_vel(_mesh.dimension() >= 2 ? coupledValue("v") : _zero),
   _w_vel(_mesh.dimension() == 3 ? coupledValue("w") : _zero),

  // Variable numberings
  _potential_var(coupled("potential")),
  _u_vel_var_number(coupled("u")),
  _v_vel_var_number(_mesh.dimension() >= 2 ? coupled("v") : libMesh::invalid_uint),
  _w_vel_var_number(_mesh.dimension() == 3 ? coupled("w") : libMesh::invalid_uint)

{}

Real
ConductionOutflowNs::computeQpResidual()
{
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  //RealVectorValue U(0, 0, 0);

  return -_test[_i][_qp]*(_grad_u[_qp]+_coefficient*_u[_qp]*(_potential_gradient[_qp]+_eqpotential_gradient[_qp]+_gradchem)
         +_coefficient*_bulkconc*std::exp(-_coefficient*_eqpotential_value[_qp])*(_potential_gradient[_qp]+_gradchem)
         -_over_diff*(_bulkconc*std::exp(-_coefficient*_eqpotential_value[_qp])+_u[_qp])*U)*_normals[_qp];
}

Real
ConductionOutflowNs::computeQpJacobian()
{
  // Derivative of the residual with respect to "u"
  RealVectorValue U1(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  //RealVectorValue U1(0, 0, 0);

  return -_test[_i][_qp]*(_grad_phi[_j][_qp]+_coefficient*_phi[_j][_qp]*(_potential_gradient[_qp]+_eqpotential_gradient[_qp]+_gradchem)
                         -_over_diff*_phi[_j][_qp]*U1)*_normals[_qp];
}

Real
ConductionOutflowNs::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Derivative of the residual with respect to "u"

  if (jvar == _potential_var)
  {
    return -_test[_i][_qp]*_coefficient*(_u[_qp]+_bulkconc*std::exp(-_coefficient*_eqpotential_value[_qp]))*_grad_phi[_j][_qp]*_normals[_qp];
  }
  else if (jvar == _u_vel_var_number)
  {
    return _over_diff*_test[_i][_qp]*(_u[_qp]+_bulkconc*std::exp(-_coefficient*_eqpotential_value[_qp]))*_phi[_j][_qp] * _normals[_qp](0);
  }
  else if (jvar == _v_vel_var_number)
  {
    return _over_diff*_test[_i][_qp]*(_u[_qp]+_bulkconc*std::exp(-_coefficient*_eqpotential_value[_qp]))*_phi[_j][_qp] * _normals[_qp](1);
  }
  else if (jvar == _w_vel_var_number)
  {
    return _over_diff*_test[_i][_qp]*(_u[_qp]+_bulkconc*std::exp(-_coefficient*_eqpotential_value[_qp]))*_phi[_j][_qp] * _normals[_qp](2);
  }
  else
    return 0.0;
}
