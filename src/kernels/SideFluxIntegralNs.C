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

#include "SideFluxIntegralNs.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<SideFluxIntegralNs>()
{
  InputParameters params = validParams<SideIntegralVariablePostprocessor>();
  params.addRequiredCoupledVar("potential", "The variable representing the pressure.");
  params.addRequiredCoupledVar("eqpotential", "The variable representing the pressure.");
  params.addRequiredParam<Real>("coefficient", "The coefficient");
  params.addRequiredParam<Real>("bulkconc", "The coefficient");
  params.addRequiredParam<RealVectorValue>("gradchem", "gradient of chemical potential");

  params.addRequiredParam<Real>("overdiff", "inverse of diffusion coefficient");

  // Coupled variables
  //params.addCoupledVar("u", "x-velocity");
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", "z-velocity"); // only required in 3D

  return params;
}

SideFluxIntegralNs::SideFluxIntegralNs(const InputParameters & parameters) :
    SideIntegralVariablePostprocessor(parameters),

  _potential_gradient(coupledGradient("potential")),
  _eqpotential_gradient(coupledGradient("eqpotential")),
  _eqpotential_value(coupledValue("eqpotential")),
  _coefficient(getParam<Real>("coefficient")),
  _bulkconc(getParam<Real>("bulkconc")),
  _gradchem(getParam<RealVectorValue>("gradchem")),

  _over_diff(getParam<Real>("overdiff")),

  // Coupled variables
  _u_vel(coupledValue("u")),
  _v_vel(_mesh.dimension() >= 2 ? coupledValue("v") : _zero),
  _w_vel(_mesh.dimension() == 3 ? coupledValue("w") : _zero)
{}

Real
SideFluxIntegralNs::computeQpIntegral()
{
  RealVectorValue U0(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  //RealVectorValue U0(0, 0, 0);

  return -(_grad_u[_qp]+_coefficient*_u[_qp]*(_potential_gradient[_qp]+_eqpotential_gradient[_qp]+_gradchem)+
  +_coefficient*_bulkconc*std::exp(-_coefficient*_eqpotential_value[_qp])*(_potential_gradient[_qp]+_gradchem)
  -_over_diff * (_bulkconc*std::exp(-_coefficient*_eqpotential_value[_qp])+_u[_qp])*U0 )*_normals[_qp];
}

