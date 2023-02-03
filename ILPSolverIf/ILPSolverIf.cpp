#include "ILPSolverIf.h"
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"

ILPSolverIf::ILPSolverIf(const SOLVER_ENUM& se) : _se(se), _t{-1}, _nvar{0}, _nrow{0}, _sol{nullptr}
{
  _solver = new OsiClpSolverInterface;
}

ILPSolverIf::~ILPSolverIf()
{
  delete static_cast<OsiClpSolverInterface*>(_solver);
  _solver = nullptr;
  delete[] _sol;
  _sol = nullptr;
}

double ILPSolverIf::getInfinity() const
{
  return static_cast<OsiClpSolverInterface*>(_solver)->getInfinity();
}

void ILPSolverIf::loadProblem(const int nvar, const int nrow, const int* starts,
        const int* indices, const double* values, const double* varlb, const double* varub,
        const double* obj, const double* rowlb, const double* rowub, const int* intvars)
{
  if (_solver) {
    if (_se == SOLVER_ENUM::Cbc) {
      _nvar = nvar;
      _nrow = nrow;
      auto sl = static_cast<OsiClpSolverInterface*>(_solver);
      sl->loadProblem(nvar, nrow, starts, indices,
          values, varlb, varub, obj, rowlb, rowub);
      for (int i = 0; i < static_cast<int>(nvar); ++i) {
        if (intvars[i]) {
          sl->setInteger(i);
        }
      }
    }
  }
}

int ILPSolverIf::solve(const int num_threads)
{
  int status{0};
  CbcModel model(*static_cast<OsiClpSolverInterface*>(_solver));
  model.setLogLevel(0);
  model.setMaximumSolutions(1000);
  model.setMaximumSavedSolutions(1000);
  if (_t > 0) model.setMaximumSeconds(_t);
  if (num_threads > 1 && CbcModel::haveMultiThreadSupport()) {
    model.setNumberThreads(num_threads);
    if (_t > 0) model.setMaximumSeconds(_t * num_threads);
    const char* argv[] = {"", "-log", "0", "-threads", std::to_string(num_threads).c_str(), "-solve"};
    status = CbcMain(6, argv, model);
  } else {
    const char* argv[] = {"", "-log", "0", "-solve"};
    status = CbcMain(4, argv, model);
  }
  status = model.secondaryStatus();
  double* var = model.bestSolution();
  if (var) {
    _sol = new double[model.getNumCols()];
    for (int i = 0; i < model.getNumCols(); ++i) _sol[i] = var[i];
  }
  return status;
}


void ILPSolverIf::writelp(char *filename, char **varnames, char **colnames)
{
  auto osiclp = static_cast<OsiClpSolverInterface*>(_solver);
  if (varnames) for (int i = 0; i < _nvar; ++i) osiclp->setColName(i, varnames[i]);
  if (colnames) for (int i = 0; i < _nrow; ++i) osiclp->setRowName(i, colnames[i]);
  osiclp->writeLp(filename);
}
