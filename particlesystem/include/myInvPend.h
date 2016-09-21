#ifndef myInvPend_h
#define myInvPend_h

#include <vector>
#include <deque>
#include <iostream>
#include <Eigen/Dense>
#include <memory>
#include "mySolver.h"
#include "myParticle.h" 
#include "myConstraint.h"
#include "myGlobConsts.h"

using namespace std;

//class that holds all refs to constraints and particles for an inverted pendulum/seaweed obj
namespace particleSystem {
	class mySolver;
	class myInvPend {
	public:
		
		myInvPend(vector<std::shared_ptr<myParticle>>& _p, int _pStIdx, int _pEndIdx, vector<std::shared_ptr<myConstraint>>& _c):
			ID(ID_gen++), p(_p),  c(_c), p_stIdx(_pStIdx), p_endIdx(_pEndIdx), partCOM(0,0,0),
			J(), Jdot(), W(), M(), q(), qdot(), Q(), feedBack(), Qhat()//, lambda()
		{
			partCOM = calcAndSetCOM();
			calcMassAndWMat();
		}

		~myInvPend() {		}//destructor	

		void buildCnstrntStruct();
		void calcMassAndWMat();
		void calcAnkleForce();
		//inline double calcAndApplyAnkleForce(int idx);                //calculate and apply appropriate ankle forces to counteract forces on constrained particle, return kp (with kd = 1/5 * deltat * kp)

		inline void calcAndApplyCnstrntFrc(bool calcCOM) {
			//calcAnkleForce();													//calc ankle forces based on derived forces

			buildCnstrntStruct();									//handle constraints here - rebuild constraint structure			
			applyConstraintForcesToSystem();									//handle collisions from contraint enforcement
			if (calcCOM) { partCOM = calcAndSetCOM(); }				//derive COM val if only 1 inv pend			 TODO find this for all inv pend in seaweed
		}

		//apply all resultant constraint forces to system (qhat)
		void applyConstraintForcesToSystem() {
			unsigned int idx = 0, qhIdx=0;
			for (unsigned int pidx = p_stIdx; pidx < p_endIdx; ++pidx) {
				p[pidx]->applyForce(Eigen::Vector3d(Qhat[qhIdx], Qhat[qhIdx + 1], Qhat[qhIdx + 2]));
				qhIdx += 3;
			}
		}//applyConstraintForcesToSystem

		Eigen::Vector3d calcAndSetCOM();

		friend ostream& operator<<(ostream& out, const myInvPend& p) {
			out << "Particle ID : " << p.ID;
			out << "\n" << endl;
			return out;
		}
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	public://variables
		static unsigned int ID_gen;
		int ID;

		Eigen::Vector3d partCOM;

	private:
		vector<std::shared_ptr<myParticle>> p;					//ref to ara of particles - only access subset defined by [p_stIdx,p_endIdx) 
		vector<std::shared_ptr<myConstraint>> c;				//array holding refs to all constraints in this seaweed

		unsigned int p_stIdx, p_endIdx;							//idx's in particle vector for the particles in this seaweed
																//constrained bead on wire matrices and vectors	n = # parts, m = num constraints
		Eigen::MatrixXd J,								//1 col per particle's x,y,z (3n) , 1 row per contraint(m)
			Jdot,							//1 col per particle's x,y,z (3n) , 1 row per contraint(m)
			W,								//diagonal inverse mass matrix, 1 3x3 diag entry per particle, rest 0 - only stores non-zero values in rows of sparse vectors (3n x 3n)
			M;								//diagonal mass matrix, 1 3x3 diag entry per particle, rest 0 - only stores non-zero values in rows of sparse vectors (3n x 3n)
		Eigen::VectorXd q,										//q vector is position of all particles
			qdot,									//velocity of all particles
			Q,										//Q is vector of forces on each particle, from force acc
			feedBack,								//complete constraint feedback term
			Qhat;									//Qhat vector is constrain forces to apply to each particle - this is what is -applied- (result of constraint calculations) to force acc
			//lambda;								//lambda is c.size() # of lagrangian multipliers, 1 per constraint


	};//class myParticle
}//namespace particleSystem
#endif//myInvPend_h