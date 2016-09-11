#include "..\include\myCollider.h"
#include <memory>

using namespace std;
namespace particleSystem{
	unsigned int myCollider::ID_gen = 0;

	//if collision type has been specified, then initialize colldier variables
	void myCollider::initCollider() {
		if (colType == FLAT) {
			buildPlaneNorm();
			findPlaneEQ();
		}
		else if (colType == SPHERE) {
			findMinMaxRadius();
		}
	}//initCollider

	//checks if particle location + deltaT partVel will cause a collision, or if particle is in contact without
	//a normal-dir component of velocity, in which case it will modify the force
	//0 vec if no collision, 1 if collision via breach, 2 if collision next timestep, 3 if need to counter force due to contact
	int myCollider::checkCollision(double deltaT, std::shared_ptr<myParticle> part) {
		Eigen::Vector3d partLoc = part->getPosition();// position[0];
		Eigen::Vector3d partVel = part->getVelocity();
		vector<Eigen::Vector3d> partVelComp;
		Eigen::Vector3d partVelNorm, partVelTan;
		Eigen::Vector3d partFAcc = part->getForceAcc();

		if (colType == FLAT) {//find distance from plane, if greater than possible displacement, then no collision, 
			double distFromBreach = (partLoc - verts[0]).dot(planeNormal);
			//else if collision, either change norm velocity to be -Krest * Velnorm, or if VelPerp == 0 then counter forces in "down" dir with equivalent force in "up" dir
			if (distFromBreach <= 0) { return 1; }						                                //immediate collision - breached
			if ((partVel.norm() != 0) && (partVel.dot(planeNormal) > 0)) { return 0; }//no collision if particle is moving in same dir as normal, so long as its position is on same side of plane as normal

			partVelComp = getPartVelNorm(partVel, planeNormal);
			partVelNorm = partVelComp[0];//planeNormal * (dot(planeNormal,partVel));
			partVelTan = partVelComp[1];//partVel - partVelNorm;
			if (distFromBreach < ((deltaT * partVel) + (.5 * partFAcc * deltaT * deltaT)).norm()) {	//means particle is within potential next-turn movement of plane
				if ((distFromBreach < epsVal) && (partVelNorm.norm() == 0)) { return 3; }				//contact with plane, no velocity in normal direction, need to counter with force in -1*planenormal dir
				return 2;																		    //potential immenent collision during next time step
			}
		}
		else if (colType == SPHERE) {//find distance
			Eigen::Vector3d partMovePoint = (partLoc + (deltaT * partVel) + (.5 * deltaT * deltaT * partFAcc));	//potential movement point for next turn of movement, to see if next turn of movement will hit wall
			if ((((partMovePoint - center).norm() < (minMaxRadius[0] + partRad)) && (!intRefl)) ||					//current location is breach 
				(((partMovePoint - center).norm() > (minMaxRadius[1] - partRad)) && (intRefl))) {
				if ((((partLoc - center).norm() < (minMaxRadius[0] + partRad)) && (!intRefl)) ||					//current location is breach 
					(((partLoc - center).norm() > (minMaxRadius[1] - partRad)) && (intRefl))) {
					return 1;
				}
				else {
					return 2;
				}
			}
			//find point on surface of sphere inline with center and partlocation
			Eigen::Vector3d sNormPartP = getSphereNormal(partLoc);				//normal through current point and center, in direction of collision surface
			Eigen::Vector3d partSpherePnt = (center + (sNormPartP * -1 * snowGlobRad));			//point on ellipsoid surface colinear with center and particle move point
			double dist2wall = (partSpherePnt - partLoc).norm();
			if (dist2wall - ((deltaT * partVel) + (.5 * deltaT * deltaT * partFAcc)).norm() > epsVal) { return 0; }
			else if (dist2wall - ((deltaT * partVel) + (.5 * deltaT * deltaT * partFAcc)).norm() > -epsVal) { return 3; }
			else { return 1; }

		}//if sphere
		return 0;
	}//checkCollision

	//calculate the normal and tangent components of velocity(or any vector) compared to a passed normal
	vector<Eigen::Vector3d> myCollider::getPartVelNorm(const Eigen::Ref<const Eigen::Vector3d>& partVel, const Eigen::Ref<const Eigen::Vector3d>& norm) {
		vector<Eigen::Vector3d> result(2, Eigen::Vector3d(0, 0, 0));
		result[0] = norm * (norm.dot(partVel));//norm dir
		result[1] = partVel - result[0];		//tan dir
		return result;
	}

	//finds minimum and maximum value of radius for ellipsoid sphere, to speed up calculations of collisions
	void myCollider::findMinMaxRadius() {
		minMaxRadius[0] = min3(radius[2], radius[0], radius[1]);
		minMaxRadius[1] = max3(radius[2], radius[0], radius[1]);
	}

	//determines the equation coefficients for a planar collider
	void myCollider::findPlaneEQ() {
		//Ax + By + Cz + D = 0
		peq[0] = planeNormal[0];		//A == norm.X
		peq[1] = planeNormal[1];		//B == norm.y
		peq[2] = planeNormal[2];		//C == norm.z
		peq[3] = (peq[0] * verts[0][0]) + (peq[1] * verts[0][1]) + (peq[2] * verts[0][2]);		//-D
	}

	//build normal of planar object
	void myCollider::buildPlaneNorm() {
		//cout<<"in build plane normal"<<endl;
		Eigen::Vector3d P0P1 = verts[1] - verts[0];
		Eigen::Vector3d P1P2 = verts[2] - verts[1];
		planeNormal = P0P1.cross(P1P2);
		planeNormal.normalize();
	}//buildNorm

	//if particle has breached planar collider somehow, move particle along appropriate normal direction until not breached anymore, and check particle's velocity(and reverse if necessary)
	Eigen::Vector3d myCollider::handlePlanarCollision(std::shared_ptr<myParticle> p, int res) {
		Eigen::Vector3d resFrc(0, 0, 0);
		Eigen::Vector3d partPos = p->getPosition();//position[0];
		double distFromBreach = (partPos - verts[0]).dot(planeNormal);
		Eigen::Vector3d partVel = p->getVelocity();
		//1 if collision via breach, 2 if collision next timestep, 3 if need to counter force due to contact
		if ((res == 2) || (partVel.dot(planeNormal) < 0)) {//going to breach next time step, or have breached and velocity is still going down - swap velocity direction
			vector<Eigen::Vector3d> partVelComp = getPartVelNorm(partVel, planeNormal);
			//partVelComp[0] *= (-1 * Krest);//reverse direction of normal velocity
			//TODO change this for verlet integrator 
			p->setVelocity( (-1 * Krest)* partVelComp[0] + partVelComp[1]);//should change dir of velocity
			// p->color = Eigen::Vector3d(1,0,0);
			if (p->solveType == VERLET) { handleVerletCol(p); }//handle reflection/velocity change by swapping old and new positions - need to scale by krest
		}//if breached and velocity going away from wall

		if (res == 1) {           //immediate breach, swap position
			//p->color = Eigen::Vector3d(0, 0, 0);
			if (distFromBreach > 0) {}//cout<<"breach error, not on wrong side of plane"<<endl;}
			else if (p->solveType == VERLET) { handleVerletCol(p); }	//handle reflection/velocity change by swapping old and new positions - need to scale by krest			
			else {//forcibly move particle to just a bit on the right side of the collider
				distFromBreach *= -(2.001);
				Eigen::Vector3d newPos((partPos + (planeNormal * (distFromBreach + epsVal))));   //reflect position up from plane by slightly more than breach amount
				//if(p->getSolveType() == GROUND){cout<<"dist from breach "<<distFromBreach<<" old position: "<<partPos<<" new position : "<<newPos<<endl;}
				p->setPosition( newPos);
				vector<Eigen::Vector3d> partAccComp = getPartVelNorm(p->getForceAcc(), planeNormal);
				Eigen::Vector3d frcTanDir = partAccComp[1].normalized();
				Eigen::Vector3d tanForce = -(1 * muFrict * (partAccComp[0].dot(planeNormal)) * frcTanDir);
				p->setForceAcc( partAccComp[0]) ;
				resFrc = -muFrict*partAccComp[1];       //apply to all particles

				p->getVelocity().setZero();// = Eigen::Vector3d(0,0,0);//partVelComp[0];//+partVelComp[1];//should change dir of velocity
			}
		}//if 1

		else if (res == 3) {          //contact
			if (p->solveType == VERLET) { handleVerletCol(p); }//handle reflection/velocity change by swapping old and new positions - need to scale by krest
		}

		if ((res == 3) || (res == 2) || (res == 1)) {                             //any contact criteria - swap normal force direction
			vector<Eigen::Vector3d> partAccComp = getPartVelNorm(p->getForceAcc(), planeNormal);
			partAccComp[0] *= -1;//reverse direction of normal acc
			p->setForceAcc(partAccComp[0] + partAccComp[1]);
			// p->applyForce((partAccComp[0] + partAccComp[1]));
		}//tangent

		return resFrc;
	}//handlePlanarBreach

	void myCollider::handleVerletCol(std::shared_ptr<myParticle> p) {
		Eigen::Vector3d tmpOldPos = p->getPosition(), tmpNewPos = p->getPosition(1);         //swapped already
		Eigen::Vector3d colPt = .5*(tmpOldPos + tmpNewPos);
		double krTmp = ((1 - Krest) * .25) + Krest;
		tmpOldPos = krTmp*(2 * (tmpOldPos - colPt.dot(planeNormal) * planeNormal) - (tmpOldPos - colPt)) + colPt;
		tmpNewPos = (2 * (tmpNewPos - colPt.dot(planeNormal) * planeNormal) - (tmpNewPos - colPt)) + colPt;
		p->setPosition(tmpNewPos);
		p->setPosition(tmpOldPos,1);
	}

	//if particle has breached planar collider somehow, move particle along appropriate normal direction until not breached anymore, and check particle's velocity(and reverse if necessary)
	void myCollider::handleSphereCollision(std::shared_ptr<myParticle> p, int res) {
		Eigen::Vector3d partPos = p->getPosition();
		Eigen::Vector3d partVel = p->getVelocity();
		Eigen::Vector3d sphereNormal = getSphereNormal(partPos);

		if (res == 2) {//if close to intersection with sphere boundary
			vector<Eigen::Vector3d> partVelComp = getPartVelNorm(partVel, sphereNormal);
			//partVelComp[0] *= (-1 * Krest);//reverse direction of normal velocity
			p->getVelocity() = ((-1 * Krest)*partVelComp[0] + partVelComp[1]);//should change dir of velocity, decrease tangent velocity for friction
		}//if about to hit collider

		else if (res == 1) {//1 if collision via breach, 2 if collision next timestep, 3 if need to counter force due to contact
			double distFromBreach = (partPos - center).norm() - (snowGlobRad - partRad);
			//cout<<"dist from breach "<<distFromBreach<<endl;
			if (((intRefl) && ((distFromBreach) < 0)) || ((!intRefl) && ((distFromBreach) > 0))) {}//cout<<"breach error, not on wrong side of sphere"<<endl;}
			else {//forcibly move particle to just a bit on the right side of the collider, reverse velocity
				distFromBreach *= (1.1);//move slightly more than breach amount
				Eigen::Vector3d newPos((partPos + (sphereNormal * (distFromBreach))));//move back into sphere
				p->setPosition(newPos);
				vector<Eigen::Vector3d> partVelComp = getPartVelNorm(partVel, sphereNormal);
				//partVelComp[0] *= -1;//reverse direction of normal velocity
				p->getVelocity() = (-1* partVelComp[0] + partVelComp[1]);//should change dir of velocity, for sphere zeroing tangent velocity
			}
		}//if 1
		else if (res == 3) //|| (res == 2) || (res == 1)) 
		{//tangent, get forceAcc and add -(forcecomponent in normal dir)
			vector<Eigen::Vector3d> partAccComp = getPartVelNorm(p->getForceAcc(), sphereNormal);
			partAccComp[0] *= -1 * Krest;//reverse direction of normal accel
			p->applyForce((partAccComp[0] + partAccComp[1]));
		}//tangent
	}//handlePlanarBreach
}//namespace particleSystem

