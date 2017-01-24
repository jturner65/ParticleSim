#ifndef myParticle_h
#define myParticle_h

#include <vector>
#include <deque>
#include <iostream>
#include <Eigen/Dense>
#include <memory>
#include "mySolver.h"
#include "myGlobConsts.h"

using namespace std;

namespace particleSystem {
	class mySolver;
	class myParticle {
		public:

			myParticle(double _mass, const Eigen::Ref<const Eigen::Vector3d>&  _pos, const Eigen::Ref<const Eigen::Vector3d>&  _velocity, const Eigen::Ref<const Eigen::Vector3d>&  _forceAcc, SolverType _solv);
			~myParticle();	

			//call only upon initiation
			inline void setOrigMass(double _m) {
				mass = _m;
				origMass = _m;
			}

			inline void reInitPartComps() {
				mass = origMass;
				color = origColor;
			}
			inline void setColor(double r, double g, double b, double a) {
				color = Eigen::Vector3d(r,g,b);
				origColor = color;
			}
			//setting these so that we can move to a single state, statedot and statedotdot vector
			inline Eigen::VectorXd getState() { return q[0]; }
			inline Eigen::VectorXd getState(int idx) { return q[idx]; }
			inline void setState(const Eigen::Ref<const Eigen::VectorXd>& _val) { q[0] << _val; }
			inline void setState(const Eigen::Ref<const Eigen::VectorXd>& _val, int idx) { q[idx] << _val; }

			inline Eigen::VectorXd getStateDot() { return qdot[0]; }
			inline Eigen::VectorXd getStateDot(int idx) { return qdot[idx]; }
			inline void setStateDot(const Eigen::Ref<const Eigen::VectorXd>& _val) { qdot[0] << _val; }
			inline void setStateDot(const Eigen::Ref<const Eigen::VectorXd>& _val, int idx) { qdot[idx] << _val; }

			inline Eigen::Vector3d getPosition() { return Eigen::Vector3d(q[0].segment<3>(0)); }
			inline Eigen::Vector3d getPosition(int idx) { return Eigen::Vector3d(q[idx].segment<3>(0)); }
			inline void modPosition(const Eigen::Ref<const Eigen::Vector3d>& _val) { q[0].segment<3>(0) += _val; }
			inline void modPosition(const Eigen::Ref<const Eigen::Vector3d>& _val, int idx) { q[idx].segment<3>(0) += _val; }
			inline void setPosition(const Eigen::Ref<const Eigen::Vector3d>& _val) { setQ(_val, 0, 0); }
			inline void setPosition(const Eigen::Ref<const Eigen::Vector3d>& _val, int idx) { setQ(_val, idx, 0); }

			inline Eigen::Vector3d getVelocity() { return Eigen::Vector3d(qdot[0].segment<3>(0));}
			inline Eigen::Vector3d getVelocity(int idx) { return Eigen::Vector3d(qdot[idx].segment<3>(0)); }
			inline void setVelocity(const Eigen::Ref<const Eigen::Vector3d>& _val) { setQDot(_val, 0, 0);  setQ(_val, 0, 3); }
			inline void setVelocity(const Eigen::Ref<const Eigen::Vector3d>& _val, int idx) { setQDot(_val, idx, 0);  setQ(_val, idx, 3); }

			inline Eigen::Vector3d getForceAcc() { return Eigen::Vector3d(qdot[0].segment<3>(3)); }
			inline Eigen::Vector3d getForceAcc(int idx) { return Eigen::Vector3d(qdot[idx].segment<3>(3)); }
			inline void setForceAcc(const Eigen::Ref<const Eigen::Vector3d>& _val) { setQDot(_val, 0, 3);  }
			inline void setForceAcc(const Eigen::Ref<const Eigen::Vector3d>& _val, int idx) { setQDot(_val, idx, 3); }

			//used to update state for 1/2 timestep calc
			inline void pushState(const Eigen::Ref<const Eigen::VectorXd>& new_st, const Eigen::Ref<const Eigen::VectorXd>& new_stDot) {q.push_front(new_st);	qdot.push_front(new_stDot);}
			//get rid of most recent temp state
			inline void popState() {		q.pop_front();		qdot.pop_front();	}

		private :
			inline void setQ(const Eigen::Ref<const Eigen::Vector3d>& _val, int idx, int segLoc) {q[idx].segment<3>(segLoc) << _val;	}
			inline void setQDot(const Eigen::Ref<const Eigen::Vector3d>& _val, int idx, int segLoc) {qdot[idx].segment<3>(segLoc) << _val;}

		public :
			//adds force to forceAcc idx 0<<
			void applyForce(const Eigen::Ref<const Eigen::Vector3d>& _force);

			void advance(const Eigen::Ref<const Eigen::VectorXd>& new_st, const Eigen::Ref<const Eigen::VectorXd>& new_stDot);

			friend ostream& operator<<(ostream& out, const myParticle& p) {
				out << "Particle ID : " << p.ID << " Mass : " << p.mass << " position : " << p.q[0].segment<3>(0) << " Velocity : " << p.qdot[0].segment<3>(0) << " Forces : " << p.qdot[0].segment<3>(3) << " Solver : " << SolverType2str[p.solveType] << " Color : " << p.color;
				out << "\n" << endl;
				return out;
			}

		private:
			deque<Eigen::VectorXd> q;
			deque<Eigen::VectorXd> qdot;
		public:
			static unsigned int ID_gen;
			int ID;
			bool inContact;
			double mass, origMass;      //handles if mass is artificially inflated
			SolverType solveType;		//{GROUND, EXP_E, MIDPOINT, RK3, RK4, IMP_E, etc }
			std::shared_ptr<mySolver> solver;			//actual  integrator used by this particle
			Eigen::Vector3d color;				//color for display of this particle
			Eigen::Vector3d origColor;				//original color for display of this particle	
			Eigen::Vector3d initPos, initVel;


	};//class myParticle
}//namespace particleSystem
#endif