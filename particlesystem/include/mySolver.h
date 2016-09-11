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

		typedef vector<Eigen::Vector3d>(mySolver::*IntegratorFuncPtr)(double, vector<Eigen::Vector3d>&, vector<Eigen::Vector3d>&);
		IntegratorFuncPtr integrator;

		mySolver(SolverType _t) :ID(++ID_gen), type(_t), lambda(2), invLam(0), lamHalf(0), hafMInvLam(0), oneMlamHalf(0) { setLambda(lambda); setIntegratorType(_t); }
		virtual ~mySolver() {}

		inline int getID() { return ID; }
		inline SolverType getType() { return type; }
		void setIntegratorType(SolverType type);

		inline double getLambda() { return lambda; }
		inline void setLambda(double _l) {
			if (_l == 0) { return; }//can't equal 0, 
			lambda = _l;  invLam = 1.0 / lambda; lamHalf = lambda/2.0, hafMInvLam = .5 - invLam, oneMlamHalf = 1.0 - lamHalf;
		}

		//for system with different solver per particle
		//vector<Eigen::Vector3d> IntegratorEvalPerPart(double deltaT, SolverType type, vector<Eigen::Vector3d>& partState, vector<Eigen::Vector3d>& partStateDot);
		//to be called directly on a per particle basis, for part 1 : _state and _stateDot are from each particle
		vector<Eigen::Vector3d> IntegrateGroundPerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot);
		vector<Eigen::Vector3d> IntegrateExp_EPerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot);
		vector<Eigen::Vector3d> IntegrateMidpointPerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot);
		vector<Eigen::Vector3d> IntegrateVerletPerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot);
		vector<Eigen::Vector3d> IntegrateRK3PerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot);
		vector<Eigen::Vector3d> IntegrateRK4PerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot);
		vector<Eigen::Vector3d> IntegrateRK4PerPartLambda(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot);
		vector<Eigen::Vector3d> IntegrateImp_EPerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot);
		vector<Eigen::Vector3d> IntegrateTrapPerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot);
		vector<Eigen::Vector3d> IntegrateNone(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot) { return _state; }

		vector<Eigen::Vector3d> Integrate(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot) {	return (this->*integrator)(deltaT, _state, _stateDot);}

		friend ostream& operator<<(ostream& out, const mySolver& s) {out << "Solver ID : " << s.ID <<  " Type : " << SolverType2str[s.type] <<std::endl;return out;}
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
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
