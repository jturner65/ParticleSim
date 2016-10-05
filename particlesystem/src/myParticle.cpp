#include "..\include\myParticle.h"

namespace particleSystem{
	unsigned int myParticle::ID_gen = 0;

	myParticle::myParticle(double _mass, const Eigen::Ref<const Eigen::Vector3d>&  _pos, const Eigen::Ref<const Eigen::Vector3d>&  _velocity, const Eigen::Ref<const Eigen::Vector3d>&  _forceAcc, SolverType _solv) :
		ID(ID_gen++), inContact(false), mass(_mass), origMass(_mass), solveType(_solv), solver(nullptr),
		color(0, 0, 0), origColor(0, 0, 0), initPos(_pos), initVel(_velocity)
		, q(2, Eigen::VectorXd(6)), qdot(2, Eigen::VectorXd(6))
	{
		//set up deques to have 2 spots
		//Eigen::VectorXd tmp(6), tmpDot(6);
		//tmp.segment<3>(0) << _pos, _velocity;
		//tmpDot.segment<3>(0) << _velocity, _forceAcc;
		solver = make_shared<mySolver>(_solv);
		q[0] << _pos, _velocity;
		q[1] << _pos, _velocity;
		qdot[0] << _velocity, _forceAcc;
		qdot[1] << _velocity, _forceAcc;
		color << (rand() / (1.0*RAND_MAX)), (rand() / (1.0*RAND_MAX)), (rand() / (1.0*RAND_MAX));
		origColor = color;

	}

	myParticle::~myParticle() {		}//destructor	

	void myParticle::applyForce(const Eigen::Ref<const Eigen::Vector3d>& _force) {
		qdot[0].segment<3>(3) += _force;
	//	forceAcc[0] += _force;
	}//applyforce 

	void myParticle::advance(const Eigen::Ref<const Eigen::VectorXd>& new_st, const Eigen::Ref<const Eigen::VectorXd>& new_stDot) {
		q.pop_back();				q.push_front(new_st);
		qdot.pop_back();			qdot.push_front(new_stDot);
	}


}//namespace partsys