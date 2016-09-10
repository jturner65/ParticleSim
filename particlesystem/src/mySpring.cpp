#include "..\include\mySpring.h"
#include "..\include\myForce.h"

namespace particleSystem{
	unsigned int mySpring::ID_gen = 0;
	//for final project for cs7492 - implicit with conjugate gradient
	void mySpring::buildSprJpJv() {
		Eigen::Vector3d delPos = (a->getPosition() - b->getPosition());
		Eigen::Matrix3d dpdpT = delPos * delPos.transpose();

		double currLen = (delPos.norm());
		if (currLen != 0) currLen = 1.0 / currLen;     //handles if currLen is 0, Jp will be 0
		dpdpT = dpdpT*(currLen*currLen);
		Jp = (dpdpT + (Id3x3 - dpdpT)*(1 - (restLen*currLen))) * Ks;
		//Jd is constant, since kd won't change
	}

	//not used, using implicit instead
	vector<Eigen::Vector3d> mySpring::calcSpringForce() {// {S_SCALAR,S_VECTOR, ATTR, SPRING};
		vector<Eigen::Vector3d> result(2, Eigen::Vector3d(0, 0, 0));
		if (restLen != 0) {		//spring with damping force
			Eigen::Vector3d aMb = (a->getPosition() - b->getPosition());
			Eigen::Vector3d lnorm = aMb.normalized();//unitlength vector of sprRestVec
			Eigen::Vector3d ldot = (a->getPosition() - a->getPosition(1));      //since using verlet, use this for velocity
			double d = ((a->getPosition(1) - b->getPosition(1)).norm());
			double KsTerm = Ks * (d - restLen);
			double KdTerm = Kd * (ldot.dot(lnorm));
			double fp = -1 * (KsTerm + KdTerm);
			result[0] = (lnorm * fp);
			result[1] = (lnorm * (-1 * fp));
		}//only add force if magnitude of distance vector is not 0
		else {
			cout << "0 restlen spring : " << ID << endl;
		}
		return result;
	}//calcForceOnParticle

}//namespace particleSystem

