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

// BCs
#include "ConductionOutflow.h"

// postprocessor
#include "SideFluxIntegralNp.h"
// Materials
#include "PackedColumn.h"

// AuxKernels
#include "Flux.h"
#include "TotalConc.h"

// functions
#include "RectFunction.h"

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

  registerAux(Flux);
  registerAux(TotalConc);
  registerBoundaryCondition(ConductionOutflow);
  registerBoundaryCondition(SideFluxIntegralNp);

  registerFunction(RectFunction);

}

// External entry point for dynamic syntax association
extern "C" void DolphinApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { DolphinApp::associateSyntax(syntax, action_factory); }
void
DolphinApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
