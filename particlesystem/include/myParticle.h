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
			ID(ID_gen++), mass(_mass), inContact(false), initPos(_pos), initVel(_velocity), solver(nullptr),
			position(1, Eigen::Vector3d(0, 0, 0)), velocity(1, Eigen::Vector3d(0, 0, 0)), forceAcc(1, Eigen::Vector3d(0, 0, 0)), 
			oldPos(1, Eigen::Vector3d(0, 0, 0)), oldVel(1, Eigen::Vector3d(0, 0, 0)), oldForceAcc(1, Eigen::Vector3d(0, 0, 0)) {
			init(_pos, _velocity, _forceAcc, _solv);
		}

		~myParticle() {}//destructor	

		void init(const Eigen::Vector3d& _pos, const Eigen::Vector3d& _velocity, const Eigen::Vector3d& _forceAcc, SolverType _solv) {
			setOrigMass(mass);
			position.push_front(_pos);
			velocity.push_front(_velocity);
			forceAcc.push_front(_forceAcc);
			oldPos.push_front(_pos);
			oldVel.push_front(_velocity);
			oldForceAcc.push_front(_forceAcc);
			solveType = _solv;
			solver = make_shared<mySolver>( _solv);

			color = Eigen::Vector3d((rand() / (1.0*RAND_MAX)), (rand() / (1.0*RAND_MAX)), (rand() / (1.0*RAND_MAX)));
			origColor = color;
		}
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

		//inline Eigen::Vector3d& getPosition() { return position[0]; }
		//inline Eigen::Vector3d& getVelocity() { return velocity[0]; }
		//inline Eigen::Vector3d& getForceAcc() { return forceAcc[0]; }

		//inline Eigen::Vector3d& getPosition(int idx) { return position[idx]; }
		//inline Eigen::Vector3d& getVelocity(int idx) { return velocity[idx]; }
		//inline Eigen::Vector3d& getForceAcc(int idx) { return forceAcc[idx]; }


		friend ostream& operator<<(ostream& out, const myParticle& p) {
			out << "Particle ID : " << p.ID << " Mass : " << p.mass << " position : "
				<< p.position[0] << " Velocity : " << p.velocity[0] << " Forces : " << p.forceAcc[0] << " Solver : " << SolverType2str[p.solveType] << " Color : " << p.color;
			out << "\n" << endl;
			return out;
		}

		//adds force to forceAcc idx 0<<
		void applyForce(const Eigen::Vector3d& _force);//{for(int idx = 0; idx<3; ++idx){forceAcc[0][idx] += _force[idx];}	cout<<"ID: "<<ID<<"applied force"<<_force<<"forceAcc[0]"<<forceAcc[0]<<endl;}
		//void applyForce(const Eigen::Vector3d& _force) { forceAcc[0] += _force; };//doesn't work in header for some reason

		inline void advance(const Eigen::Vector3d& _p, const Eigen::Vector3d& _v, const Eigen::Vector3d& _a) {
			oldPos.pop_back();          oldPos.push_front(position.front());            //save old pos
			position.pop_back();		position.push_front(_p);	//add new first element, move the rest of the elements along 1 step

			oldVel.pop_back();          oldVel.push_front(velocity.front());             //save old vel
			velocity.pop_back();    	velocity.push_front(_v);

			oldForceAcc.pop_back();     oldForceAcc.push_front(forceAcc.front());
			forceAcc.pop_back();        forceAcc.push_front(_a);
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

		deque<Eigen::Vector3d> oldPos;		//previous time step position is idx currIDX : is a deque to hold previous terms
		deque<Eigen::Vector3d> oldVel;		//previous time step velocity is idx currIDX
		deque<Eigen::Vector3d> oldForceAcc;   //previous time step forceAcc val is idx currIDX
		Eigen::Vector3d initPos, initVel;

		deque<Eigen::Vector3d> position;		//current position is idx currIDX : is a deque to hold previous terms
		deque<Eigen::Vector3d> velocity;		//current velocity is idx currIDX
		deque<Eigen::Vector3d> forceAcc;		//current forceAcc val is idx currIDX

	private:
		deque<Eigen::VectorXd> q;
		deque<Eigen::VectorXd> qdot;


	};//class myParticle
}//namespace particleSystem
#endif