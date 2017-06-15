#include "DolphinApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

// Kernels
#include "Potential.h"
#include "EqPotential.h"
#include "Concentration.h"
#include "ConcTimeDerivative.h"

// ns kernels
#include "ConvectionVel.h"
#include "INSMomentumForceLaplaceForm.h"
#include "INSMomentumForceLaplaceFormRZ.h"
#include "INSMomentumForceTractionForm.h"
#include "INSMomentumForceTractionFormRZ.h"

#include "InterfaceDiffusion.h"

// BCs
#include "ConductionOutflow.h"
#include "ConductionOutflowNs.h"

// postprocessor
#include "SideFluxIntegralNp.h"
#include "SideFluxIntegralNs.h"

// Materials
#include "PackedColumn.h"

// AuxKernels
#include "Flux.h"
#include "FluxNs.h"
#include "TotalConc.h"

// functions
#include "RectFunction.h"
#include "RectFunctionExp.h"
#include "RectFunctionExp2.h"

template<>
InputParameters validParams<DolphinApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

DolphinApp::DolphinApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  DolphinApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  DolphinApp::associateSyntax(_syntax, _action_factory);
}

DolphinApp::~DolphinApp()
{
}

// External entry point for dynamic application loading
extern "C" void DolphinApp__registerApps() { DolphinApp::registerApps(); }
void
DolphinApp::registerApps()
{
  registerApp(DolphinApp);
}

// External entry point for dynamic object registration
extern "C" void DolphinApp__registerObjects(Factory & factory) { DolphinApp::registerObjects(factory); }
void
DolphinApp::registerObjects(Factory & factory)
{
  registerKernel(Potential);
  registerKernel(EqPotential);
  registerKernel(Concentration);
  registerKernel(ConcTimeDerivative);

  registerKernel(INSMomentumForceLaplaceForm);
  registerKernel(INSMomentumForceLaplaceFormRZ);

  registerKernel(INSMomentumForceTractionForm);
  registerKernel(INSMomentumForceTractionFormRZ);

  registerKernel(ConvectionVel);

  registerAux(Flux);
  registerAux(FluxNs);

  registerAux(TotalConc);

  // Interface kernels
  registerInterfaceKernel(InterfaceDiffusion);

  registerBoundaryCondition(ConductionOutflow);
  registerBoundaryCondition(ConductionOutflowNs);

  registerBoundaryCondition(SideFluxIntegralNp);
  registerBoundaryCondition(SideFluxIntegralNs);

  registerFunction(RectFunction);
  registerFunction(RectFunctionExp);
  registerFunction(RectFunctionExp2);

}

// External entry point for dynamic syntax association
extern "C" void DolphinApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { DolphinApp::associateSyntax(syntax, action_factory); }
void
DolphinApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
