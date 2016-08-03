#include "..\include\myParticle.h"

namespace particleSystem{
	unsigned int myParticle::ID_gen = 0;

	void myParticle::applyForce(const Eigen::Vector3d& _force) {
		forceAcc[0] += _force;
		//oldForceAcc[1] = Eigen::Vector3d(forceAcc[0]);
	}//applyforce


}//namespace partsys