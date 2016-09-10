#include "..\include\myParticle.h"

namespace particleSystem{
	unsigned int myParticle::ID_gen = 0;

	void myParticle::applyForce(const Eigen::Vector3d& _force) {
		qdot[0].segment<3>(6) += _force;
		forceAcc[0] += _force;
	}//applyforce 


}//namespace partsys