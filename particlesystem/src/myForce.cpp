#include "..\include\myForce.h"

namespace particleSystem{
	unsigned int myForce::ID_gen = 0;

	//old version, leave in place
	vector<Eigen::Vector3d> myForce::calcForceOnParticle(std::shared_ptr<myParticle> _p1, std::shared_ptr<myParticle> _p2, double d, std::shared_ptr<myForce> force) {// {S_SCALAR,S_VECTOR, ATTR, SPRING};
		vector<Eigen::Vector3d> result(2, Eigen::Vector3d(0, 0, 0));
		ForceType fTyp = force->ftype;
		Eigen::Vector3d l, v_l;
		switch (fTyp) {
			case S_SCALAR:	{
				result[0] = (force->constVec * _p1->mass);
				break; }
			case S_VECTOR:	{
				result[0] = (_p1->velocity[0] * force->constVal);//vel 
				break; }
			case ATTR: {//attractor, uses two particles, 1st constant, 
				l = _p1->position[0] - _p2->position[0];
				double lmag = l.norm();
				if (lmag > epsVal) {		//attractor force - constVal (negative) * m1 * m2 * lnorm/ lmag*lmag
					double m1 = _p1->mass, m2 = _p2->mass;
					Eigen::Vector3d lnorm = l.normalized();//unitlength vector of l
					double fp = -1 * force->constVal * m1 * m2 / (lmag * lmag);
					result[0] = (lnorm * fp);
					result[1] = (lnorm * (-1 * fp));
				}//only add force if magnitude of distance vector is not 0
				break; }
			case REPL: {//repulsive force, uses two particles, 1st constant, opposite sign as attractor 
				l = _p1->position[0] - _p2->position[0];
				double lmag = l.norm();
				//if(lmag < 1){		//repulsive force -> -( - constVal * m1 * m2 * lnorm/ lmag*lmag) - only kick into play when closer than 1
				double m1 = _p1->mass, m2 = _p2->mass;
				Eigen::Vector3d lnorm = l.normalized();//unitlength vector of l
				double fp = -1 * force->constVal * m1 * m2 / (lmag * lmag);
				result[0] = (-1 * (lnorm * fp));
				result[1] = (-1 * (lnorm * (-1 * fp)));
				//}//only add force if magnitude of distance vector is not 0
				break; }
			case DAMPSPRING:{//damped spring - not sure if going to use, but what the hey - dependent on old length (need ldot vector)
				l = _p1->position[0] - _p2->position[0];
				double lmag = l.norm();
				if (lmag != 0) {		//spring with damping force
					Eigen::Vector3d lnorm = l.normalized();//unitlength vector of l
					Eigen::Vector3d lprime = l - (_p1->position[1] - _p2->position[1]);		//lprime - time derivative, subtract old length vector from new length vector ?
					double KsTerm = force->constVal * (lmag - d);
					double KdTerm = force->constVal2 * (lprime.dot(l));
					double fp = -1 * (KsTerm + KdTerm);
					result[0] = (lnorm * fp);
					result[1] = (lnorm * (-1 * fp));
				}//only add force if magnitude of distance vector is not 0
				break; }
			case DSPR_THETABAR:{//damped spring to represent ankle
				l = _p1->position[0] - _p1->initPos;
				v_l = _p1->velocity[0] - _p1->initVel;
				double lmag = l.norm();
				if (lmag != 0) {		//spring with damping force
					Eigen::Vector3d lnorm = l.normalized();//unitlength vector of l
					//Eigen::Vector3d lprime = l - (_p1->position[1] - _p2->position[1]);		//lprime - time derivative, subtract old length vector from new length vector ?
					double KpTerm = force->constVal * (lmag - d);
					double KdTerm = force->constVal2 * (v_l.dot(l));
					double fp = -1 * (KpTerm + KdTerm);
					result[0] = (lnorm * fp);
					result[1] = (lnorm * (-1 * fp));
				}//only add force if magnitude of distance vector is not 0
				break; }
			default: {	break; }

		}//switch
		return result;
	}//calcForceOnParticle

}//namespace particleSystem

