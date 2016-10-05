#ifndef mySolver_h
#define mySolver_h

#include <vector>
#include <string>
#include <iostream>
#include <Eigen/Dense>
#include "myGlobConsts.h"
#include "myParticle.h"
#include <memory>

using namespace std;

namespace particleSystem{
	class mySolver {
	public:

		typedef Eigen::VectorXd(mySolver::*IntegStateFuncPtr)(double, const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&);
		IntegStateFuncPtr integrator;

		mySolver(SolverType _t);
		virtual ~mySolver();

		//mySolver(SolverType _t) :ID(++ID_gen), type(_t), lambda(2), invLam(0), lamHalf(0), hafMInvLam(0), oneMlamHalf(0) { setLambda(lambda); setIntegratorType(_t);}
		//virtual ~mySolver() {}

		inline int getID() { return ID; }
		inline SolverType getType() { return type; }

		inline double getLambda() { return lambda; }
		inline void setLambda(double _l) {
			if (_l == 0) { ++_l; }//can't equal 0, 
			lambda = _l;  invLam = 1.0 / lambda; lamHalf = lambda/2.0, hafMInvLam = .5 - invLam, oneMlamHalf = 1.0 - lamHalf;
		}

		void setIntegratorType(SolverType type);
		Eigen::VectorXd IntegrateGroundPerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot);
		Eigen::VectorXd IntegrateExp_EPerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot);
		Eigen::VectorXd IntegrateMidpointPerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot);
		Eigen::VectorXd IntegrateVerletPerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot);
		Eigen::VectorXd IntegrateRK3PerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot);
		Eigen::VectorXd IntegrateRK4PerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot);
		Eigen::VectorXd IntegrateRK4PerPartLambda(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot);
		Eigen::VectorXd IntegrateImp_EPerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot);
		Eigen::VectorXd IntegrateTrapPerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot);
		Eigen::VectorXd IntegrateNone(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot) { return _state; }


		Eigen::VectorXd Integrate(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot) { return (this->*integrator)(deltaT, _state, _stateDot); }


		friend ostream& operator<<(ostream& out, const mySolver& s) {out << "Solver ID : " << s.ID <<  " Type : " << SolverType2str[s.type] <<std::endl;return out;}
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	protected:
		Eigen::VectorXd IntegrateExp_EFrcUpd(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot);
	public:
		static unsigned int ID_gen;
		int ID;
		SolverType type;
		double oldDelT;                         //for verlet integration in case time step changes
	private:
		double lambda, invLam, lamHalf, hafMInvLam, oneMlamHalf;                          //value for rk4 general form
	
	};//class mySolver
}//particleSystem namespace
#endif
