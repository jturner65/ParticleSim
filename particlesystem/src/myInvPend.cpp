#include "..\include\myInvPend.h"

namespace particleSystem {
	unsigned int myInvPend::ID_gen = 0;

	//used for inv pendulum
	void myInvPend::calcAnkleForce() {
		//torque to force on particle = tau @ "ankle" / (len of leg * sin theta)
		//theta == angle from vertical = 
		//tau_net = I * thetDotDot

		//lagrangian formulation
		//p1 is base of constraint, p2 is particle to receive force 

		for (unsigned int cIdx = 0; cIdx < c.size(); ++cIdx) {
			double kp = 0, cLen = c[cIdx]->c_Dist;
			//apply correcting force to particle 2 idx - applied at particle 1 (treating like "ankle" joint)
			int pIdx = c[cIdx]->p2Idx;
			Eigen::Vector3d
				//ankleDisp = (c[cIdx]->getCurAnklePoint() - c[cIdx]->getInitAnklePoint()),				//displacement of ankle (base particle/anchor point) from init position
				//curInitDisp = p[pIdx]->initPos + ankleDisp,												//current displacement of initial/target position of pendulum body
				dispFromInit = c[cIdx]->getPendDistFromInit(),
				tmpThet = dispFromInit +
				(.1 * c[cIdx]->getCurPendVel()); //vector of positions and velocities
			double lenFacc = p[pIdx]->getForceAcc().norm(),
				lenThet = tmpThet.norm()/cLen;
			kp = (0 == lenThet ? 0 : lenFacc * 1.1 / lenThet);
			//kp *= 1;
			Eigen::Vector3d newForce = (0 == kp ? Eigen::Vector3d(0, 0, 0) : -kp * lenThet * (p[pIdx]->getForceAcc().normalized()));
			p[pIdx]->applyForce(newForce);       //should be 0
		}
	}//calcAndApplyAnkleForce
	 //find COM of particle configuration
	void myInvPend::calcAndSetCOM() {
		partCOM.setZero();
		//int numParts = p.size();
		double totMass = 0;
		for (unsigned int pidx = p_stIdx; pidx < p_endIdx; ++pidx) {
			partCOM += p[pidx]->getPosition() * p[pidx]->mass;
			totMass += p[pidx]->mass;
		}
		partCOM = partCOM / totMass;
	}//calcAndSetCOM

	//precalculate M (mass) matrix and W matrix - only do when configuration changes
	void myInvPend::calcMassAndWMat() {
		int num3Parts = (p_endIdx - p_stIdx) * 3;
		W = Eigen::MatrixXd(num3Parts, num3Parts);
		M = Eigen::MatrixXd(num3Parts, num3Parts);
		W.setIdentity();
		M.setIdentity();
		int nIdx = 0;
		for (unsigned int idx = p_stIdx; idx < p_endIdx; ++idx) {
			double val1 = 1.0 / p[idx]->mass;
			for (int idxOff = 0; idxOff < 3; ++idxOff) {
				W(nIdx + idxOff, nIdx + idxOff) = val1;
				M(nIdx + idxOff, nIdx + idxOff) = p[idx]->mass;//1/( 1.0 * val1);
			}
			nIdx += 3;
		}
	}//calcMassAndWMat


	void myInvPend::buildCnstrntStruct() {
		//init J, Jdot matrices; lambda, C vecs, Cdot vecs
		unsigned int numParts = (p_endIdx - p_stIdx) * 3;
		unsigned int numCnstrnts = c.size();

		J.setZero(numCnstrnts, numParts);
		Jdot.setZero(numCnstrnts, numParts);
		q.setZero(numParts);
		qdot.setZero(numParts);
		Q.setZero(numParts);
		feedBack.setZero(numCnstrnts);

		int nIdx = 0;							//index into part-related vals is nIdx, constraint-related vals is mIdx
												//first set up q, qdot, qdotdot vectors
		for (unsigned int idx = p_stIdx; idx < p_endIdx; ++idx) {
			//cout << "calc constraint for idx : " << idx << " part id : " << p[idx]->ID << " qIDX : " << (nIdx) << endl;
			q.segment<3>(nIdx) = p[idx]->getPosition();
			qdot.segment<3>(nIdx) = p[idx]->getVelocity();
			Q.segment<3>(nIdx) = p[idx]->getForceAcc();
			nIdx += 3;		//increment by 3 for each particle
		}
		int p1nIDX, p2nIDX;
		//next build constraint values for each constraint endpoint
		for (unsigned int cIdx = 0; cIdx < numCnstrnts; ++cIdx) {
			p1nIDX = 3 * (c[cIdx]->p1Idx - p_stIdx);
			p2nIDX = 3 * (c[cIdx]->p2Idx - p_stIdx);

			feedBack(cIdx) = c[cIdx]->ks * c[cIdx]->calcCVal(p[c[cIdx]->p1Idx], p[c[cIdx]->p2Idx]) + c[cIdx]->kd * c[cIdx]->calcCDotVal(p[c[cIdx]->p1Idx], p[c[cIdx]->p2Idx]);
			//col for each particle is 3*constraint's pidx, row for each particle is constraint idx 
			J.row(cIdx).segment<3>(p1nIDX) = c[cIdx]->calcPartialCP1(p[c[cIdx]->p1Idx], p[c[cIdx]->p2Idx]);
			Jdot.row(cIdx).segment<3>(p1nIDX) = c[cIdx]->calcPartialCdotP1(p[c[cIdx]->p1Idx], p[c[cIdx]->p2Idx]);
			J.row(cIdx).segment<3>(p2nIDX) = c[cIdx]->calcPartialCP1(p[c[cIdx]->p2Idx], p[c[cIdx]->p1Idx]);
			Jdot.row(cIdx).segment<3>(p2nIDX) = c[cIdx]->calcPartialCdotP1(p[c[cIdx]->p2Idx], p[c[cIdx]->p1Idx]);
		}

		Eigen::MatrixXd JW = J*W;
		Eigen::MatrixXd Jtrans = J.transpose();
		//lambda = (JW * Jtrans).ldlt().solve(-Jdot*qdot - JW * Q - feedBack);
		Qhat = Jtrans * (JW * Jtrans).ldlt().solve(-Jdot*qdot - JW * Q - feedBack);			//constraint force applied

	}//void buildConstraintStructure();


}//namespace partsys