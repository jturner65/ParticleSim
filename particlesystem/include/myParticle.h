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
		
		myParticle(double _mass, const Eigen::Vector3d&  _pos, const Eigen::Vector3d&  _velocity, const Eigen::Vector3d&  _forceAcc, SolverType _solv) :
			ID(ID_gen++), inContact(false), mass(_mass), origMass(_mass), solveType(_solv),solver(make_shared<mySolver>(_solv)), 
			color((rand() / (1.0*RAND_MAX)), (rand() / (1.0*RAND_MAX)), (rand() / (1.0*RAND_MAX))), origColor(0,0,0),initPos(_pos), initVel(_velocity),
			position(2,_pos), velocity(2,_velocity), forceAcc(2,_forceAcc)
			,q(2, Eigen::VectorXd(6)),qdot(2, Eigen::VectorXd(6))
		{
			//init(_pos, _velocity, _forceAcc);

			//solveType = _solv;
			//solver = make_shared<mySolver>(_solv);

			//color = Eigen::Vector3d((rand() / (1.0*RAND_MAX)), (rand() / (1.0*RAND_MAX)), (rand() / (1.0*RAND_MAX)));
			origColor = color;

		}

		~myParticle() {}//destructor	

		//call only upon initiation
		void setOrigMass(double _m) {
			mass = _m;
			origMass = _m;
		}

		void reInitPartComps() {
			mass = origMass;
			//color = origColor;
		}
		void setColor(double r, double g, double b, double a) {
			color = Eigen::Vector3d(r,g,b);
			origColor = color;
		}
		//setting these so that we can move to a single state, statedot and statedotdot vector
		inline Eigen::VectorXd& getState() { return q[0]; }
		inline Eigen::VectorXd& getState(int idx) { return q[idx]; }
		inline void setState(const Eigen::VectorXd& _val) { q[0] << _val; }
		inline void setState(const Eigen::VectorXd& _val, int idx) { q[idx] << _val; }

		inline Eigen::VectorXd& getStateDot() { return qdot[0]; }
		inline Eigen::VectorXd& getStateDot(int idx) { return qdot[idx]; }
		inline void setStateDot(const Eigen::VectorXd& _val) { qdot[0] << _val; }
		inline void setStateDot(const Eigen::VectorXd& _val, int idx) { qdot[idx] << _val; }

		//inline Eigen::Vector3d& getPosition() { return Eigen::Vector3d(q[0].segment<3>(0)); }
		//inline Eigen::Vector3d& getPosition(int idx) { return Eigen::Vector3d(q[idx].segment<3>(0)); }
		//inline void setPosition(const Eigen::Vector3d& _val) { q[0].segment<3>(0) << _val; }
		//inline void setPosition(const Eigen::Vector3d& _val, int idx) { q[idx].segment<3>(0) << _val; }

		inline Eigen::Vector3d& getPosition() { return position[0]; }
		inline Eigen::Vector3d& getPosition(int idx) { return position[idx]; }
		inline void setPosition(const Eigen::Vector3d& _val) { position[0] << _val; }
		inline void setPosition(const Eigen::Vector3d& _val, int idx) { position[idx] << _val; }

		inline Eigen::Vector3d& getVelocity() { return velocity[0]; }
		inline Eigen::Vector3d& getVelocity(int idx) { return velocity[idx]; }
		inline void setVelocity(const Eigen::Vector3d& _val) { velocity[0] << _val; qdot[0].segment<3>(0) << _val; q[0].segment<3>(3) << _val; }
		inline void setVelocity(const Eigen::Vector3d& _val, int idx) { velocity[idx] << _val; qdot[idx].segment<3>(0) << _val; q[idx].segment<3>(3) << _val;}

		inline Eigen::Vector3d& getForceAcc() { return forceAcc[0]; }
		inline Eigen::Vector3d& getForceAcc(int idx) { return forceAcc[idx]; }
		inline void setForceAcc(const Eigen::Vector3d& _val) { forceAcc[0] << _val;qdot[0].segment<3>(3) << _val;  }
		inline void setForceAcc(const Eigen::Vector3d& _val, int idx) { forceAcc[idx] << _val; qdot[idx].segment<3>(3) << _val; }


		friend ostream& operator<<(ostream& out, const myParticle& p) {
			out << "Particle ID : " << p.ID << " Mass : " << p.mass << " position : "
				<< p.position[0] << " Velocity : " << p.velocity[0] << " Forces : " << p.forceAcc[0] << " Solver : " << SolverType2str[p.solveType] << " Color : " << p.color;
			out << "\n" << endl;
			return out;
		}

		//adds force to forceAcc idx 0<<
		void applyForce(const Eigen::Vector3d& _force);//{for(int idx = 0; idx<3; ++idx){forceAcc[0][idx] += _force[idx];}	cout<<"ID: "<<ID<<"applied force"<<_force<<"forceAcc[0]"<<forceAcc[0]<<endl;}

		inline void advance(Eigen::VectorXd& new_st,Eigen::VectorXd& new_stDot) {
			//q.pop_back();				q.push_front(new_st);
			//qdot.pop_back();			qdot.push_front(new_stDot);

			position.pop_back();		position.push_front(new_st.segment<3>(0));	//add new first element, move the rest of the elements along 1 step
			velocity.pop_back();    	velocity.push_front(new_st.segment<3>(3));
			forceAcc.pop_back();        forceAcc.push_front(new_stDot.segment<3>(3));
		}


		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
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

	private:
		deque<Eigen::Vector3d> position;		//current position is idx currIDX : is a deque to hold previous terms
		deque<Eigen::Vector3d> velocity;		//current velocity is idx currIDX
		deque<Eigen::Vector3d> forceAcc;		//current forceAcc val is idx currIDX

		deque<Eigen::VectorXd> q;
		deque<Eigen::VectorXd> qdot;

		//void init(const Eigen::Vector3d& _pos, const Eigen::Vector3d& _velocity, const Eigen::Vector3d& _forceAcc) {
		//	//set up deques to have 2 spots
		//	//Eigen::VectorXd tmp(6), tmpDot(6);
		//	//tmp.segment<3>(0) << _pos;
		//	//tmp.segment<3>(3) << _velocity;
		//	//tmpDot.segment<3>(0) << _velocity;
		//	//tmpDot.segment<3>(3) << _forceAcc;
		//	//q[0] = (tmp);
		//	//q[1] = (tmp);
		//	//qdot[0] = (tmpDot);
		//	//qdot[1] = (tmpDot);

		//	position.push_front(Eigen::Vector3d(_pos));
		//	position.push_front(Eigen::Vector3d(_pos));
		//	velocity.push_front(Eigen::Vector3d(_velocity));
		//	velocity.push_front(Eigen::Vector3d(_velocity));
		//	forceAcc.push_front(Eigen::Vector3d(_forceAcc));
		//	forceAcc.push_front(Eigen::Vector3d(_forceAcc));
		//	//oldPos.push_front(_pos);
		//	//oldVel.push_front(_velocity);
		//	//oldForceAcc.push_front(_forceAcc);
		//}


	};//class myParticle
}//namespace particleSystem
#endif