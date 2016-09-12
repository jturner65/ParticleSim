#include "..\include\mySystem.h"
#include "..\include\MyParticleWorld.h"
#include <vector>

namespace particleSystem{
	unsigned int mySystem::ID_gen = 0;
	//constraint spring constant
	double mySystem::globKp = 1000;

    void mySystem::applyForceVecToAllPartsMsSpr(const Eigen::Ref<const Eigen::Vector3d>& frc){
        int numParts = p.size();
		for(int i = 0; i<numParts-1; ++i){
            p[i]->applyForce(frc);
		}//for every particle
	}//apply forces to all particles
	//handle particle-particle collisions no repel force
	void mySystem::handlePartPartCollision() {
		unsigned int psize = p.size();
		for (unsigned int pIdx = 0; pIdx < psize - 1; ++pIdx) {
			for (unsigned int qIdx = pIdx + 1; qIdx < psize; ++qIdx) {
				double dist = (p[pIdx]->getPosition() - p[qIdx]->getPosition()).squaredNorm();
				if (dist < partSqColDist) {// check if particles within potential collision dist of each other
					if (checkParticlePairForCollision(p[pIdx], p[qIdx], dist)) {					//check if headed in same direction
						//if so, swap velocities
						Eigen::Vector3d pVel(p[pIdx]->getVelocity());
						p[pIdx]->setVelocity(p[qIdx]->getVelocity());
						p[qIdx]->setVelocity(pVel);
					}//if collision, then swap velocities
				}
			}//qIdx = this is to pick both particles to compare - need to set both particles' values in the compare -  this should only make any pairing one time		
		}//don't need to handle last particle from p idx - will get handled with every other check via q idx		
	}//handlePartPartCollision

	//handle particle-particle collisions
	void mySystem::handlePartPartRepelCollision() {//adds repeling force
		myForce::ID_gen = f.size();
		//std::shared_ptr<myForce> tmpForceRepl = std::make_shared<myForce>("Repulsive Dummy",2);			//dummy force
		Eigen::Vector3d tmp(0, 0, 0);
		vector<Eigen::Vector3d> tmpRes(2,Eigen::Vector3d(0,0,0));
		unsigned int psize = p.size();
		for(unsigned int pIdx = 0; pIdx < psize -1; ++pIdx){
			for(unsigned int qIdx = pIdx + 1; qIdx < psize; ++qIdx){
				tmp = (p[pIdx]->getPosition() - p[qIdx]->getPosition());
				double dist = tmp.squaredNorm();
				//calcReplForceOnParticles(std::shared_ptr<myParticle> _p1, std::shared_ptr<myParticle> _p2, const Eigen::Ref<const Eigen::Vector3d>& diffVec, double kVal)
				tmpRes = myForce::calcReplForceOnParticles(p[pIdx], p[qIdx], tmp, 2);		//apply repulsion force
				p[pIdx]->applyForce(tmpRes[0]);
				p[qIdx]->applyForce(tmpRes[1]);
				if(dist < partSqColDist){// check if particles within potential collision dist of each other
					if (checkParticlePairForCollision(p[pIdx], p[qIdx], dist)){					//check if headed in same direction
							//if so, swap velocities
						Eigen::Vector3d pVel(p[pIdx]->getVelocity());
						p[pIdx]->setVelocity(p[qIdx]->getVelocity());
						p[qIdx]->setVelocity(pVel);
					}//if collision, then swap velocities
				}
			}//qIdx = this is to pick both particles to compare - need to set both particles' values in the compare -  this should only make any pairing one time		
		}//don't need to handle last particle from p idx - will get handled with every other check via q idx		
	}//handlePartPartCollision

	//checks if two particles are going to collide with one another in the next time step, if so, swap the components of their velocities in the directions they collide
	bool mySystem::checkParticlePairForCollision(std::shared_ptr<myParticle> partA, std::shared_ptr<myParticle> partB, double dist){
		//make one particle stationary
		Eigen::Vector3d combVel = partA->getVelocity() - partB->getVelocity();				//part b is "stationary", and part a has all velocity
		Eigen::Vector3d velDelT = combVel * deltaT;
		Eigen::Vector3d partDiffs = partA->getPosition() - partB->getPosition();
		Eigen::Vector3d partBtoADir = (partA->getPosition() - partB->getPosition()).normalized();	//vector from particle b to particle a - if velocity vector not close to this direction, then no collision
		if ((velDelT.squaredNorm() < (dist * dist)) || (partBtoADir.dot(combVel) >= 0)) { return false; }	//particles not going fast enough in the appropriate direction to get closer to each other - moving apart (sqr calc faster than sqrt)
																									//or a's velocity in wrong direction => perp or in the same general direction as vector from b to a
		double colRad = 2 * partRad;											//minimum distance within which a collision will occur - treat like the radius of a circle that partA ray is going to possibly intersect
		//Eigen::Vector3d partAPos = partA->position[0], partBPos = partB->position[0], partDiffs = partAPos - partBPos;
		//determine ray-sphere intersect, where ray is partA.pos + (combVelDir * t) and sphere is ctred at partB.pos and radius colRad - if t works out to be < delta t then collision occurs
		double A = combVel.squaredNorm();	//Vx^2 + Vy^2 + Vz^2
		double B = 2 * (combVel.dot(partDiffs));		//partB pos = sphere center, partA pos, ray origin
		double C = (partDiffs.squaredNorm()) - (colRad*colRad);
		//quadratic eq - discr first
		double discr = (B*B - (4*A*C));
		if (discr < 0) {return false;}		//imaginary roots - no intersection
		double t1 = ((-1*B) + sqrt(discr))/(2.0 * A);
		double t2 = ((-1*B) - sqrt(discr))/(2.0 * A);
		if((t1 >= 0) && (t1 < deltaT) || (t2 >= 0) && (t2 < deltaT)) {return true;}			//intersection happens
		return false;				//no intersection
	}//checkParticlePairForCollision

	//check all collisions - particle to particle and square/spherical boundary(ec)
	bool mySystem::handlePartCldrCollision(){//check for each collider for all particles - check returns 0 - no collision possible, 1- breach, 2- collision due to velocity, 3- collision due to contact force
		int res, psize = p.size(), cldrsize = colliders.size();
        Eigen::Vector3d resFrc(0,0,0), tmpFrc(0,0,0);
	
		for (int i = 0; i < cldrsize; ++i) {
			for (int partIDX = 0; partIDX < psize; ++partIDX) {
				res = colliders[i]->checkCollision(deltaT, p[partIDX]);
				if(res > 0){					
					if(colliders[i]->colType == FLAT){				
                        tmpFrc = colliders[i]->handlePlanarCollision(p[partIDX], res);			
                        resFrc += tmpFrc;
                    }
					else if(colliders[i]->colType == SPHERE){		colliders[i]->handleSphereCollision(p[partIDX], res);			}
					else { cout<<" unknown collision type : "<<colliders[i]->colType<<endl;}
					//if(p[partIDX]->getID() % 5 == 0){cout<<"res of collision : "<<res<<"\n"<<*this<<endl;}
				}
			}//for each particle
		}//for each collider		
		if (spr.size() != 0) { applyForceVecToAllPartsMsSpr(resFrc); }//this handles friction force for being on the ground - need to find alternate mechanism for "rolling" mass spring TODO
        if(resFrc.norm() >0){return true;}
        return false;
 	}//check for collisions
    
    //initialize/reinitialize when all particles/constraints/colliders have been set, or when any have been added or changed
	void mySystem::initMassSystem() {
		if (flags[useMassMat]) {		calcMassAndWMat();		}
	}

	//void mySystem::handleGlobeTimeStep(double fMult, bool& dragged, const Eigen::Ref<const Eigen::Vector3d>& msdrgVal0, const Eigen::Ref<const Eigen::Vector3d>& msdrgVal1) {
	//	//fluidBox->resetOldVals();
	//	if ((dragged) || ((shakeVal).norm() > 0.001)) {
	//		addShakeForceToFluid(fMult, msdrgVal0, msdrgVal1); 		dragged = false;			//setting false keeps force from compounding
	//	}
	//	fluidBox->myFluidBoxTimeStep();					//timestep snowglobe	
	//	//fluidBox->resetOldVals();
	//	calcFluidForceForAllParticles();				//calculate fluid forces
	//	shakeVal = Eigen::Vector3d(0, 0, 0);
	//}

    //precalculate M (mass) matrix and W matrix - only do when configuration changes
    void mySystem::calcMassAndWMat(){
        int num3Parts = 3*p.size();
        W = Eigen::MatrixXd(num3Parts,num3Parts);
        M = Eigen::MatrixXd(num3Parts,num3Parts);
		W.setIdentity();
		M.setIdentity();
        int nIdx = 0;
        for (unsigned int idx = 0; idx < p.size(); ++idx){
            double val1 = 1.0/p[idx]->mass;
			for(int idxOff = 0; idxOff < 3; ++idxOff){
				W(nIdx+idxOff,nIdx+idxOff) = val1;
				M(nIdx+idxOff,nIdx+idxOff) = p[idx]->mass;//1/( 1.0 * val1);
			}
            nIdx +=3;
        }
    }//calcMassAndWMat
    
    //use implicit method to solve for spring system
    void mySystem::solveImpEuler4MassSpring(){
        //calculate every spring's jacobians
        int numSprings = spr.size(), numParts = p.size();
        for(int i =0; i<numSprings; ++i){         spr[i]->buildSprJpJv();      }        //resets Jp for new position info

        //building system to solve A x  = b to find using conjugate gradient method(maybe - A might be easily inverted), where 
        //x = dv, the linearized acceleration at t+1 needed for implicit calc of v_t+1 - once we have this, we can find pos_t+1 using pos_t+1 = pos_t + deltaT * v_t+1
        //A = M - dt*Jv - dt^2*Jp, a 3nx3n matrix 
        //b = dt*(forceVec + dt*Jp*v0), a 3n vector       
        int numParts3 = 3*numParts;                
        Eigen::VectorXd b(numParts3), p0(numParts3), v0(numParts3), f0(numParts3);
		b.setZero(); p0.setZero(); v0.setZero(); f0.setZero();
        Eigen::MatrixXd A(numParts3, numParts3);
		A.setZero();
        A = M;      //mass matrix
        int sIdx = 0, pIdx = 0, vIdx = 0;
		for (int i = 0; i<numParts; ++i) {     //get all positions, velocities, forces for every particle
			p0.segment<3>(sIdx) = p[i]->getPosition();
			v0.segment<3>(sIdx) = p[i]->getVelocity();
			f0.segment<3>(sIdx) = p[i]->getForceAcc();
			sIdx += 3;
		}
        int aIdx, bIdx;
        double delT2 = deltaT*deltaT;
        Eigen::Vector3d tmp1;
        Eigen::Matrix3d mtmp2, mtmp3;
        //build A matrix : M - deltaT*Jv - deltaT^2 * Jp
		Eigen::VectorXd delTJv(numParts3);    delTJv.setZero();
		Eigen::VectorXd delT2Jp(numParts3);    delT2Jp.setZero();
		Eigen::VectorXd dfdxV0(numParts3);    dfdxV0.setZero();

        for(int i=0; i<numSprings; i++ ){//mult df/dx by rel velocity between two parts
            tmp1 = spr[i]->Jp*(spr[i]->a->getVelocity() - spr[i]->b->getVelocity());
            mtmp2 = deltaT * spr[i]->Jv;
            mtmp3 = delT2 * spr[i]->Jp;
            aIdx = 3 * spr[i]->a->ID;
            bIdx = 3 * spr[i]->b->ID;

            dfdxV0[aIdx] -= tmp1[0];          dfdxV0[aIdx+1] -= tmp1[1];          dfdxV0[aIdx+2] -= tmp1[2];
            dfdxV0[bIdx] += tmp1[0];          dfdxV0[bIdx+1] += tmp1[1];          dfdxV0[bIdx+2] += tmp1[2];
            //build A matrix
            for(int j = 0; j < 3; ++j){
                for(int k = 0; k<3; ++k){
                    A(aIdx+k,bIdx+j) -= mtmp2(k,j);
                    A(aIdx+k,bIdx+j) -= mtmp3(k,j);
                    A(bIdx+k,aIdx+j) -= mtmp2(k,j);
                    A(bIdx+k,aIdx+j) -= mtmp3(k,j);
                }
            }
        }
        b = deltaT * (f0 + (deltaT*dfdxV0));
        Eigen::VectorXd dvNew = calcConjGrad(b, A, f0);            //TODO solve this

        Eigen::VectorXd v1 = v0 + deltaT * dvNew;
		Eigen::VectorXd p1 = p0 + deltaT * v1;
		Eigen::VectorXd nstate (6), nstateDot(6);
        pIdx = 0;        vIdx = 0;
		for (int i = 0; i< numParts; ++i) {           //set results from forward integration
			//newPos[0] = p1[pIdx++];          newPos[1] = p1[pIdx++];            newPos[2] = p1[pIdx++];
			//newVel[0] = v1[vIdx++];          newVel[1] = v1[vIdx++];            newVel[2] = v1[vIdx++];
			nstate << p1.segment<3>(pIdx), v1.segment<3>(vIdx);
			nstateDot << v1.segment<3>(vIdx), Eigen::Vector3d(0, 0, 0);
			p[i]->advance(nstate, nstateDot);
			pIdx += 3; vIdx += 3;
		}

        //for(int i = 0; i< numParts; ++i){           //set results from forward integration
        //    newPos[0] = p1[pIdx++];          newPos[1] = p1[pIdx++];            newPos[2] = p1[pIdx++];
        //    newVel[0] = v1[vIdx++];          newVel[1] = v1[vIdx++];            newVel[2] = v1[vIdx++]; 
        //    p[i]->advance(newPos, newVel, Eigen::Vector3d(0,0,0));
        //}
    }//implicit solver

	Eigen::VectorXd mySystem::calcConjGrad(const Eigen::Ref<const Eigen::VectorXd>& b, const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::VectorXd>& f0) {
        int numParts = 3*p.size(), iters=0;
		Eigen::VectorXd v(numParts); v.setZero();
        Eigen::VectorXd res = b - A*v;
        Eigen::VectorXd d = res;
        Eigen::VectorXd q = b;
        double epsNew = res.dot(res), epsOld, alpha, beta, eps = .000001;
        //iterative loop
        while((iters < CONJGRAD_MAXITERS) && (epsNew > eps)){
            q = A*d;
            alpha  = epsNew/(d.dot(q));
            v = v + (alpha * q);
            res = res - alpha*q;
            epsOld = epsNew;
			epsNew = res.dot(res);
            beta = epsNew/epsOld;
            d = res + beta*d;
            iters++;
        }
        return v;
    }

    //move the central particle sufficiently to see if we can get this thing to move
    void mySystem::dispCenterPart(Eigen::Vector3d& dir){
        dir[1] = 0;
        if(dir.norm()<.0001){return;}
        //for(int i =0; i<numSprings; ++i){        spr[i]->calcNewRLenInDir(dir);     }        //resets Jp for new position info
        int numParts = p.size(), numSprings = spr.size();
		Eigen::Vector3d ctrOrigPos = calcAndSetCOM(0, p.size() - 1), dispVec;
        std::shared_ptr<myParticle> ctr = p[numParts-1];
        Eigen::Vector3d newPos = ctrOrigPos + .01*dir;
        
        //find particle closest to this one
        double minDist = 99999, d;
        int minIdx = -1;
        for(int i=0; i<numParts-2; ++i){ 
            d = (newPos - p[i]->getPosition()).norm();
            if(d<minDist){             minDist = d;             minIdx  = i;          }
        }
        dispVec = p[minIdx]->getPosition() - ctrOrigPos;
        Eigen::Vector3d nDir = dir.normalized();
		//ctr->position[0] = ctrOrigPos + (.95*dispVec.norm()) * nDir;//95% toward the closest particle from central spot
		ctr->setPosition(ctrOrigPos + (.95*dispVec.norm()) * nDir);//95% toward the closest particle from central spot
		p[minIdx]->mass = ctr->mass;
		initMassSystem();
        for(int i =0; i<numSprings; ++i){        spr[i]-> reCalcRestLen();      }   //recalculate rest length

    }
    void mySystem::resetCtrPointMassSpr(){
        int numParts = p.size(), numSprings = spr.size();
		Eigen::Vector3d ctrOrigPos = calcAndSetCOM(0, p.size() - 1), dispVec;
        p[numParts-1]->setPosition( ctrOrigPos);
        for(int i =0; i<numSprings; ++i){        spr[i]-> reCalcRestLen();      }   //recalculate rest length
    }
	//find COM of particle configuration
	Eigen::Vector3d mySystem::calcAndSetCOM(unsigned int stIdx, unsigned int endIdx) {
		Eigen::Vector3d tmpCOM(0, 0, 0);
		//int numParts = p.size();
		double totMass = 0;
		for (unsigned int pidx = stIdx; pidx < endIdx; ++pidx) {
			tmpCOM += p[pidx]->getPosition() * p[pidx]->mass;
			totMass += p[pidx]->mass;
		}
		tmpCOM = tmpCOM / totMass;
		return tmpCOM;
	}//calcAndSetCOM

    //satisfy all constraints via positional modification of particles-  use verlet integration to handle fwd int, and not have velocity issues
    void mySystem::satisfyPosConstraints(){
		unsigned int numCnstrnts = c.size(), numIters = numPosCnstIters;
        for(unsigned int i =0; i < numIters; ++i){
            for (unsigned int cIdx = 0; cIdx < numCnstrnts; ++cIdx){
                c[cIdx]->satisfyPosConstraint();
            }        
        }
    }//satisfyPosConstraints

	//this will apply all single-particle forces to all particles, clearing out the force accumulator first
	void mySystem::applyForcesToSystem() {
		vector<Eigen::Vector3d> result;
		//clearParticleForceAcc();		//clear all force acc's - cleared after deriv eval
		int numForces = f.size(),
			numParts = p.size();
		for (int fidx = 0; fidx < numForces; ++fidx) {
			for (int pidx = 0; pidx < numParts; ++pidx) {
				//second particle and d ignored for single particle forces
				result = myForce::calcForceOnParticle(p[pidx], p[pidx], 0, f[fidx]);//do not use d for single particle forces
				p[pidx]->applyForce(result[0]);
			}//for each force
		}//for all particles
	}//applyForcesToSystem

	//apply spring forces to all particles
	void mySystem::applySpringForcesToSystem() {
		int numSprings = spr.size();
		vector<Eigen::Vector3d> result;
		for (int i = 0; i<numSprings; ++i) {
			result = spr[i]->calcSpringForce();
			spr[i]->a->applyForce(result[0]);
			spr[i]->b->applyForce(result[1]);
		}
	}

	//x,y,z are edge bounds of flat collider at specific y location
	void mySystem::buildGndCollider(double kRest, double muFrict, double x, double y, double z, std::string& name, const Eigen::Ref<const Eigen::Vector3d>& drawLoc){//, const Eigen::Ref<const Eigen::Vector3d>& gndLoc) {
		vector<Eigen::Vector3d> tmpVec(4, Eigen::Vector3d(0, 0, 0));
		tmpVec[0] = Eigen::Vector3d(-x, y, 0);
		tmpVec[1] = Eigen::Vector3d(x, y, 0);
		tmpVec[2] = Eigen::Vector3d(0, y, -z);
		tmpVec[3] = Eigen::Vector3d(0, y, z);
		myCollider::ID_gen = colliders.size();
		colliders.push_back(std::make_shared<myCollider>(name, drawLoc, tmpVec, false));
		colliders.back()->Krest = kRest;
		colliders.back()->muFrict = muFrict;        //friction force
	//	ground = gndLoc;
	}

	void mySystem::buildGlobeCollider(double krest, double muFrict, double rad, double distFromGlb) {
		myCollider::ID_gen = colliders.size();
		colliders.push_back(std::make_shared<myCollider>("SnowGlobePart2", Eigen::Vector3d(0, 0, -(rad + distFromGlb)),
			Eigen::Vector3d(0, 0, -(rad + distFromSnowGlobe)),
			Eigen::Vector3d(rad, rad, rad), true));
		colliders.back()->Krest = krest;
		colliders.back()->muFrict = muFrict;
	}

	void mySystem::_buildBaseDefForces(double dCoeff) {
		std::string str1(name);
		str1.append("_Force_Grav");
		myForce::ID_gen = f.size();
		f.push_back(std::make_shared<myForce>(str1, gravVec, S_SCALAR));
		if (0 != dCoeff) {	
			std::string str2(name);
			str2.append("_Fluid_Drag");
			f.push_back(std::make_shared<myForce>(str2, dCoeff, S_VECTOR));
		}
	}

	void mySystem::buildAndSetCnstrnts(int part1IDX, int part2IDX, double rad, Eigen::Vector3d& center) {
		string name;
		if (rad == -1) { rad = dist2Parts(part1IDX, part2IDX);	name = "CircleBarbellConstraint"; }			//part1IDX != part2IDX;
		else {				name = "CircleConstraint";	}											//part1IDX == part2IDX;
		myConstraint::ID_gen = c.size();
		c.push_back(std::make_shared<myConstraint>(name, rad, mySystem::globKp, mySystem::globKp * 1.5 * deltaT, C_Circular, center));
		c.back()->setP1((p[part1IDX]), part1IDX);
		c.back()->setP2((p[part2IDX]), part2IDX);
		c.back()->useAnchor = part1IDX == part2IDX;
	}

	//build single inverted pendulum
	void mySystem::buildInvPend(Eigen::Vector3d& sLoc) {
		myParticle::ID_gen = p.size();
		for (int pIdx = 0; pIdx < 4; ++pIdx) {//add 4 new particles for single "inv pendulum"

			addParticle(1.0, Eigen::Vector3d(sLoc[0], sLoc[1] + pIdx, sLoc[2]), RK4);
			p.back()->mass = 1 + (5 - pIdx / 10.0);
			if (0 == pIdx) {
				buildAndSetCnstrnts(p.size() - 1, p.size() - 1, 1, Eigen::Vector3d(sLoc[0], sLoc[1] - 1, sLoc[2]));	            //first recent particle
				c.back()->drawCnstrPath = false;                                               //don't draw circle for first constraint of inverted pendulum
			}
			else {
				buildAndSetCnstrnts(p.size() - 2, p.size() - 1, -1, p[p.size() - 1]->getPosition());
			}
		}
	}//buildInvPend

	void mySystem::buildRollerCoasterConstraints(int id, double rad) {
		buildAndSetCnstrnts(id, id, rad, Eigen::Vector3d(2, 2, -5));
		buildAndSetCnstrnts(id, id, rad, Eigen::Vector3d(-2, 2, -5));
		buildAndSetCnstrnts(id, id, rad, Eigen::Vector3d(2, -2, -5));
		buildAndSetCnstrnts(id, id, rad, Eigen::Vector3d(-2, -2, -5));
	}

	//calculate the effects of the fluid sim on the particles in the globe - to take the place of rigged forces in shaking routine
	void mySystem::calcFluidForceForAllParticles() {
		Eigen::Vector3d fluidFrcVal(0, 0, 0);
		double shakeValMag, logShakeVal;
		for (unsigned int pidx = 0; pidx < p.size(); ++pidx) {
			//Eigen::Vector3d pPos = p[pidx]->position[0];
			fluidFrcVal = fluidBox->getVelAtCell(p[pidx]->getPosition());
			shakeValMag = fluidFrcVal.norm();
			if (shakeValMag > .00001) {
				logShakeVal = log(shakeValMag) / log(1.1);
				//Eigen::Vector3d applyShakeVal = (norm(tmpResult)) * 2 * log(shakeValMag);//(logShakeVal < 6 ? 6 + log(1 + logShakeVal) : 1+log(logShakeVal)) ;
				Eigen::Vector3d applyShakeVal = (fluidFrcVal.normalized()) * ((1 + logShakeVal) > 6 ? 6 + log(1 + logShakeVal) : (1 + logShakeVal));
				p[pidx]->applyForce(applyShakeVal);// * (colliders[0].getDistFromCenter(p[pidx]->position[0], 1)/(1.0 * (snowGlobRad))));
				//addDraggingForce(applyShakeVal, colliders[0], 1);	//collider 0 is the snow globe sphere
			}
		}//for all particles
		shakeVal = Eigen::Vector3d(0, 0, 0);
	}//calcFluidForceForAllParticles

	//add external force to extremal cells of fluid - walls tangent to plane which force is applied to
	//apply force in this plane to walls, dotted by their normals
	void mySystem::addShakeForceToFluid(double fMult, const Eigen::Ref<const Eigen::Vector3d>& msdrgVal0, const Eigen::Ref<const Eigen::Vector3d>& msdrgVal1) {
		shakeVal = fMult*(msdrgVal1 - msdrgVal0);
		int sX = fluidBox->numCellX, sY = fluidBox->numCellY, sZ = fluidBox->numCellZ;

		double  xCoord = (msdrgVal0[0] > 350 ? sX - 1 : (msdrgVal0[0] < -350 ? 0 : (sX * (msdrgVal0[0] + 350) / 700))),			//change values to be center of display
			yCoord = (msdrgVal0[1] > 350 ? sY - 1 : (msdrgVal0[1] < -350 ? 0 : (sY * (msdrgVal0[1] + 350) / 700)));

		//cout << "fluid shake force val : " << evec3dToStr(shakeVal) << "mouse0 : " << evec3dToStr(msdrgVal0) << " x : " << xCoord << " y : " << yCoord << endl;
		Eigen::Vector3d tmpVec(0, 0, 0);
		//for( int yidx = 0; yidx < sY; ++yidx){
		double scZ = 1.0 / sZ;
		for (int zidx = 0; zidx < sZ / 2; ++zidx) {
			//tmpVec = Eigen::Vector3d(shakeVal) * scZ;
			tmpVec = Eigen::Vector3d(shakeVal * zidx / sZ);
			fluidBox->myFluidBoxAddForce(Eigen::Vector3d(xCoord, yCoord, zidx), tmpVec);
		}
		for (int zidx = sZ / 2; zidx < sZ; ++zidx) {
			//tmpVec = Eigen::Vector3d(shakeVal) * scZ;
			tmpVec = Eigen::Vector3d(shakeVal * (sZ - zidx) / sZ);
			fluidBox->myFluidBoxAddForce(Eigen::Vector3d(xCoord, yCoord, zidx), tmpVec);
		}
	}//addForceToFluid		

	void mySystem::addControlMS2(double& tsCntr, int numMSMsclSpr) {
		// 0:{0,7}, 1:{0,8}, 2:{2,7}, 3:{2,8}, 4:{4,7}, 5:{4,8}, //left side - idx 4&5 are "up" springs
		//6:{1,11}, 7:{1,12}, 8:{3,11}, 9:{3,12}, 10:{4,11}, 11:{4,12}     //right side springs idx 10&11 are "up" springs
		//lift front legs as back legs lower
		//pull in lifted legs as low legs push out
		double mult = .15;
		double driveFront = mult*sin(tsCntr), driveBack = mult*sin(tsCntr - (.125 * PI)), driveFast = mult*sin(2 * tsCntr - (.333 * PI));
		for (int i = 0; i<4; ++i) {
			spr[i]->setRestLenByMult(driveBack);      //back legs foot control
			spr[i + 6]->setRestLenByMult(driveFront);   //front legs foot control
		}

		for (int i = 4; i<6; ++i) {
			spr[i]->setRestLenByMult(-driveFront);      //back legs foot control
			spr[i + 6]->setRestLenByMult(-driveBack);   //front legs foot control
		}
		for (int i = 12; i < numMSMsclSpr; ++i) {
			spr[i]->setRestLenByMult(0);
		}
		tsCntr += deltaT;
	}

	Eigen::Vector3d mySystem::calcDragForce(Eigen::Vector3d& sValP) {
		double shakeValMag = (sValP).norm() + 1;
		double logShakeVal = log(shakeValMag);
		Eigen::Vector3d dirShakeVal = ((shakeValMag > 0) ? ((sValP).normalized()) : Eigen::Vector3d(0, 0, 0));
		//cout << "svalp:" << evec3dToStr(sValP) << " mag :" << shakeValMag << " log : " << logShakeVal << " Dir : " << evec3dToStr(dirShakeVal) << endl;
		Eigen::Vector3d applyShakeVal = dirShakeVal * 2 * logShakeVal;//(logShakeVal < 6 ? 6 + log(1 + logShakeVal) : 1+log(logShakeVal)) ;
		return applyShakeVal;
	}//calcAndAddForce

	//handle adding forces to tinkertoys
	void mySystem::addForcesToTinkerToys(double fMult, bool msDragged, const Eigen::Ref<const Eigen::Vector3d>& msDragDif) {
		if (msDragged) { shakeVal = fMult*(msDragDif); }
		if (shakeVal.norm() > 0) {
			Eigen::Vector3d res = calcDragForce(shakeVal);
			//cout << "res : " << evec3dToStr(res) << endl;
			addDraggingForceTT(res, 1);	//collider 0 is the snow globe sphere	
			shakeVal *= .99;
		}
	}//addForcesToTinkerToys

	//add force to center particle in mass-spring ball
	void mySystem::addForcesToMassSpringCtrl(double fMult, bool msDragged, const Eigen::Ref<const Eigen::Vector3d>& msDragDif) {
		if (msDragged) { shakeVal = fMult*(msDragDif); }
		if (shakeVal.norm() > 0) {
			Eigen::Vector3d res = calcDragForce(shakeVal);
			//addDraggingForceTT(res, 1);	//collider 0 is the snow globe sphere
			//p[numRTMsclSpr]->applyForce(res);
			dispCenterPart(res);
			shakeVal *= .95;
		}
	}//addForcesToTinkerToys

	//add force from dragging on tinkertoys
	void mySystem::addDraggingForceTT(Eigen::Vector3d& force, double mult) {
		for (unsigned int pidx = 0; pidx < p.size(); ++pidx) {
			p[pidx]->applyForce(force);
		}//for all particles
	}//addDraggingForceTT

	//used for inv pendulum
	void mySystem::calcAnkleForce() {
		kpAra = vector<double>(p.size());
		for (unsigned int pidx = 0; pidx < p.size(); ++pidx) {
			kpAra[pidx] = calcAndApplyAnkleForce(pidx);
		}
	}//calcAndApplyAnkleForce

	//for inv pend calculate and apply appropriate ankle forces to counteract forces on constrained particle, return kp (with kd = 1/5 * deltat * kp)
	double mySystem::calcAndApplyAnkleForce(int pidx) {
		//assume constraints length 1
		double kp = 0;          
		Eigen::Vector3d tmpThet = (p[pidx]->initPos - p[pidx]->getPosition() + (.1 * p[pidx]->getVelocity())); //vector of positions and velocities
		double lenFacc = p[pidx]->getForceAcc().norm(),
			lenThet = tmpThet.norm();
		kp = (0 == lenThet ? 0 : lenFacc * 1.1 / lenThet);
		//kp *= 1;
		Eigen::Vector3d newForce = (0 == kp ? Eigen::Vector3d(0, 0, 0) : -kp * lenThet * (p[pidx]->getForceAcc().normalized()));
		p[pidx]->applyForce(newForce);       //should be 0
		return kp;
	}//calcAndApplyAnkleForce

	void mySystem::invokeSolverDerivEval() {
		Eigen::VectorXd state(12), stateDot(12), res(6);
		state.setZero();
		stateDot.setZero();
		Eigen::VectorXd tmp(6), tmp2(6);
		unsigned int numParts = p.size();
		//int solverTypeToUse = (usePartSolver ? p[idx]->solveType : solveType);
		for (unsigned int idx = 0; idx < numParts; ++idx) {
			state.segment<6>(0) << p[idx]->getState();
			state.segment<6>(6) << p[idx]->getState(1);

			stateDot.segment<6>(0) << p[idx]->getStateDot();
			stateDot.segment<6>(6) << p[idx]->getStateDot(1);

			stateDot.segment<3>(3) /= p[idx]->mass;
			stateDot.segment<3>(9) /= p[idx]->mass;

			res = p[idx]->solver->Integrate(deltaT, state, stateDot);
			
			tmp << res.segment<6>(0);
			tmp2 << res.segment<3>(3), Eigen::Vector3d(0, 0, 0);//force isn't integrated
			p[idx]->advance(tmp, tmp2);        //clears forces from force acc
		}//for each particle
	}//invokeSolverDerivEval
	
	bool mySystem::handlePauseDrawClick(int& clickOnPartIDX, const Eigen::Ref<const Eigen::Vector3d>& mClickLoc, bool idxIs5) {
		bool sceneModded = false;
		bool clickOnPart = false;
		bool releasedOnPart = false;
		unsigned int partidx, cnstidx;
		bool clickOnCnstrnt = false;
		partidx = 0;
		while ((!clickOnPart) && (partidx < p.size())) {
			//cout << "part to check click loc : " << evec3dToStr(p[partidx]->position[0]) << endl;
			if ((mClickLoc - p[partidx]->getPosition()).norm() <= 2 * partRad) { clickOnPart = true;	sceneModded = true; break; }
			++partidx;
		}//while
		cnstidx = 0;
		//check existing path constraints - if clicked on existing constraint and not on existing particle, create particle here and create constraint for it
		while ((!clickOnCnstrnt) && (cnstidx < c.size())) {
			double dist = abs((mClickLoc - c[cnstidx]->anchorPoint).norm() - c[cnstidx]->c_Dist);
			if ((c[cnstidx]->useAnchor) && (dist <= 2 * partRad)) {
				clickOnCnstrnt = true;
				sceneModded = true;
				break;
			}
			++cnstidx;
		}//while
		//cout << "click on part : " << clickOnPart << endl;
		//if clicked on existing constraint and not on existing particle, create particle here and put it on this constraint by creating a new one for it right here
		if ((!clickOnPart) && (clickOnCnstrnt)) {
			myParticle::ID_gen = p.size();
			addParticle(1.0, Eigen::Vector3d(mClickLoc), RK4);
			int idx = p.size() - 1;
			Eigen::Vector3d cCtr(c[cnstidx]->anchorPoint);
			buildAndSetCnstrnts(idx, idx, c[cnstidx]->c_Dist, cCtr);//put last particle as member of new concentric constraint
			buildCnstrntStruct(idxIs5);
			clickOnPartIDX = -1;
			sceneModded = true;
			flags[buildDragCnstrnt] = false;		//building constraint from a particle to another particle or space (path constraint) - do not make a path constraint if particle is already on a path constraint						
		}
		else if (!clickOnPart) {//means we're not on a particle and we're not on an existing constraint - create a new particle here
			clickOnPartIDX = -1;
			flags[buildDragCnstrnt] = false;		//building constraint from a particle to another particle or space (path constraint) - do not make a path constraint if particle is already on a path constraint
			myParticle::ID_gen = p.size();
			addParticle(1.0, Eigen::Vector3d(mClickLoc), RK4);
			sceneModded = true;
		}
		else {				//means we have clicked on an existing particle so we need to set up building and drawing constriant link between this particle and either space (new path constraint) or another particle(tinkertoy constraint)
			clickOnPartIDX = partidx;
			flags[buildDragCnstrnt] = true;		//building constraint from a particle to another particle or space (path constraint) - do not make a path constraint if particle is already on a path constraint
			if (clickOnCnstrnt) { flags[canBuildDragPathCnstrnt] = false; }
			else { flags[canBuildDragPathCnstrnt] = true; }
		}
		return sceneModded;
	}

	bool mySystem::handlePauseDrawRel(int& clickOnPartIDX, const Eigen::Ref<const Eigen::Vector3d>& mRelLoc, bool idxIs5) {
		bool sceneModded = false;
		bool clickOnPart = false;
		bool releasedOnPart = false;
		unsigned int partidx;// , cnstidx;
		bool clickOnCnstrnt = false;
		partidx = 0;

		//just dragged the mouse to a location, either in a particle, or in space.  
		//if in space, check to make sure initial click was not on a path-constrained particle, and if so then build path constraint with radius = drag distance
		//if on a particle, construct a constraint to this particle from intitial source particle
		if (flags[buildDragCnstrnt]) {
			while ((!releasedOnPart) && (partidx < p.size())) {
				//cout << "part to check rel loc : " << evec3dToStr(p[partidx]->position[0]) << endl;
				if ((mRelLoc - p[partidx]->getPosition()).norm() <= 2 * partRad) { releasedOnPart = true; break; }
				++partidx;
			}//while
			//cout << "release : released on particle : " << releasedOnPart << " part idx : " << partidx << endl;
			if (!releasedOnPart) {//means we're not on a particle when we released - create path constraint origin here, if able to, with particle clickOnPartIDX
				if ((flags[canBuildDragPathCnstrnt]) && (clickOnPartIDX != -1)) {
					Eigen::Vector3d cCtr(mRelLoc);
					buildAndSetCnstrnts(clickOnPartIDX, clickOnPartIDX, (p[clickOnPartIDX]->getPosition() - mRelLoc).norm(), cCtr);//put last particle as member of new concentric constraint
					buildCnstrntStruct(idxIs5);
					flags[buildDragCnstrnt] = false;		//building constraint from a particle to another particle or space (path constraint) - do not make a path constraint if particle is already on a path constraint
				}
			}
			else {//means we are on a particle when we released, build constraint here from clickOnPartIDX to this new particle
				if (clickOnPartIDX != -1) {//if actually clicked on a particle the first time
					Eigen::Vector3d cCtr(mRelLoc);
					buildAndSetCnstrnts(clickOnPartIDX, partidx, -1, cCtr);//put last particle as member of new concentric constraint
					buildCnstrntStruct(idxIs5);
					flags[buildDragCnstrnt] = false;		//building constraint from a particle to another particle or space (path constraint) - do not make a path constraint if particle is already on a path constraint
				}
			}
			sceneModded = true;     //ending a drag constraint
		}//if build drag constraint
		flags[buildDragCnstrnt] = false;
		flags[canBuildDragPathCnstrnt] = false;
		clickOnPartIDX = -1;
		return sceneModded;
	}

	//solve matrix eqs, once they are built, for lambda
	void mySystem::calcConstraintForces(){
		Eigen::MatrixXd JW = J*W;
		Eigen::MatrixXd Jtrans = J.transpose();

		lambda = (JW * Jtrans).ldlt().solve(-Jdot*qdot - JW * Q - CVal - CDotVal);

		Qhat = Jtrans * lambda;				//constraint force applied

	}//void calcConstraintForces();

	void mySystem::buildCnstrntStructJumper() {
		//init J, Jdot matrices; lambda, C vecs, Cdot vecs
		int numParts = p.size() * 3;
		int numCnstrnts = c.size();
		Eigen::MatrixXd tmpJ(numCnstrnts, numParts);			//1 col per particle's x,y,z (3n) , 1 row per constraint(m)
		Eigen::MatrixXd tmpJdot(numCnstrnts, numParts);			//1 col per particle's x,y,z (3n) , 1 row per constraint(m)
		tmpJ.setZero();
		tmpJdot.setZero();
		Eigen::VectorXd tmpCVal(numCnstrnts);							//c evaluated
		Eigen::VectorXd tmpCDotVal(numCnstrnts);						//c dot evaluated
		tmpCVal.setZero();
		tmpCDotVal.setZero();

		Eigen::VectorXd tmpq(numParts);							//q vector is position of all particles
		Eigen::VectorXd tmpqdot(numParts);						//qdot vector is velocity of all particles
		Eigen::VectorXd tmpqdotdot(numParts);							//Q vector is forces on each particle - qdoubledot will be WQ (matrix mult)
		tmpq.setZero();
		tmpqdot.setZero();
		tmpqdotdot.setZero();

		Eigen::Vector3d
			tmpVec4, tmpVec5,
			partLocNext;

		int nIdx = 0;							//index into part-related vals is nIdx, constraint-related vals is mIdx
		for (unsigned int idx = 0; idx < p.size(); ++idx) {

			int mIdx = 0;
			double cnstDist = 99999999999;
			int cnstIDX = -1;
			partLocNext = p[idx]->getPosition() + (deltaT *  p[idx]->getVelocity()) + (.5 * deltaT * deltaT * p[idx]->getForceAcc());//approx location of particle next timestep if unconstrained
																																	 //for jumping constraints
			double tmpDist;
			for (unsigned int cIdx = 0; cIdx < numCnstrnts; ++cIdx) {
				//if((p[idx]->ID == c[cIdx]->p1ID) || ((p[idx]->ID != c[cIdx]->p1ID) && (p[idx]->ID ==  c[cIdx]->p2ID))){		//if particle is p1 or p2 in constraint cnstr
				if ((p[idx]->ID == c[cIdx]->p1ID) && (p[idx]->ID == c[cIdx]->p2ID)) {		//if particle is p1 or p2 in constraint cnstr
					tmpDist = (partLocNext - c[cIdx]->anchorPoint).norm();
					if (cnstDist > tmpDist) {//find closest constraint to where unconstrained particle will move to next turn - use this constraint in calculation below
						cnstDist = tmpDist;
						cnstIDX = cIdx;
					}
				}
			}//for each constraint find constraint that is closest to this particle
			 //cout<<"part : "<<idx<<" : cnstrnt : "<<cnstIDX<<" size of constrt struct : "<<c.size()<<" next loc : "<<partLocNext<<endl;
			if (cnstIDX == -1) {
				continue;
			}//means particle has no constraint, continue with particle loop

			tmpq.segment<3>(nIdx) = p[idx]->getPosition();
			tmpqdot.segment<3>(nIdx) = p[idx]->getVelocity();
			tmpqdotdot.segment<3>(nIdx) = p[idx]->getForceAcc();

			//populate per constraint/per particle structures

			for (unsigned int cIdx = 0; cIdx < numCnstrnts; ++cIdx) {
				if (cnstIDX == cIdx) { //closest path constraint
					if ((p[idx]->ID == c[cIdx]->p1ID)) {		//if particle is p1 in constraint cnstr
																//cout<<"\np1 calc:"<<endl;
						tmpVec4 = c[cIdx]->calcPartialCP1(p[idx], p[c[cIdx]->p2Idx]);	//path constraint
						tmpVec5 = c[cIdx]->calcPartialCdotP1(p[idx], p[c[cIdx]->p2Idx]);
					}
					else if ((p[idx]->ID != c[cIdx]->p1ID) && (p[idx]->ID == c[cIdx]->p2ID)) {		//if particle is p1 in constraint cnstr
																									//cout<<"\np2 calc:"<<endl;
						tmpVec4 = c[cIdx]->calcPartialCP2(p[c[cIdx]->p1Idx], p[idx]);
						tmpVec5 = c[cIdx]->calcPartialCdotP2(p[c[cIdx]->p1Idx], p[idx]);
					}
				}
				else if ((p[idx]->ID != c[cIdx]->p2ID) && (p[idx]->ID == c[cIdx]->p1ID)) {		//if particle is p1 in constraint cnstr and not a path constraint
																								//cout<<"\np2 calc:"<<endl;
					tmpVec4 = c[cIdx]->calcPartialCP1(p[idx], p[c[cIdx]->p2Idx]);
					tmpVec5 = c[cIdx]->calcPartialCdotP1(p[idx], p[c[cIdx]->p2Idx]);
				}
				else if ((p[idx]->ID != c[cIdx]->p1ID) && (p[idx]->ID == c[cIdx]->p2ID)) {		//if particle is p2 in constraint cnstr and not a path constraint
																								//cout<<"\np2 calc:"<<endl;
					tmpVec4 = c[cIdx]->calcPartialCP2(p[c[cIdx]->p1Idx], p[idx]);
					tmpVec5 = c[cIdx]->calcPartialCdotP2(p[c[cIdx]->p1Idx], p[idx]);
				}
				else {
					//cout<<"\nno calc"<<endl;
					tmpVec4.setZero();
					tmpVec5.setZero();
				}
				tmpCVal[cIdx] = c[cIdx]->ks * c[cIdx]->calcCVal(p[c[cIdx]->p1Idx], p[c[cIdx]->p2Idx]);
				tmpCDotVal[cIdx] = c[cIdx]->kd * c[cIdx]->calcCDotVal(p[c[cIdx]->p1Idx], p[c[cIdx]->p2Idx]);
				//SPD IMPLEMENTATION TODO
				//tmpCVal[cIdx] = c[cIdx].ks * c[cIdx].calcCValSPD(p[c[cIdx].p1Idx],p[c[cIdx].p2Idx], deltaT);		//.8 = ks
				//tmpCDotVal[cIdx] = c[cIdx].kd * c[cIdx].calcCDotValSPD(p[c[cIdx].p1Idx],p[c[cIdx].p2Idx], deltaT);//.8 = kd

				tmpJ.row(cIdx).segment<3>(nIdx) = tmpVec4;
				tmpJdot.row(cIdx).segment<3>(nIdx) = tmpVec5;
			}
			nIdx += 3;		//increment by 3 for each particle
		}//for each particle

		 //now need to remove from J matrix every row that is completely 0, otherwise JWJt will not be invertible, and need to remove equivalent row from Jdot, component from lambda vector and component from cval and cdotval vectors
		int matLen = 0;		//number of rows in new matrix
		Eigen::MatrixXd tmpJ2(numCnstrnts, numParts), tmpJDot2(numCnstrnts, numParts);			//will then shrink these once we have all 0-impact constraints removed
		Eigen::VectorXd tmpCVal2(numCnstrnts), tmpCDotVal2(numCnstrnts);
		Eigen::VectorXd tmpJRow, tmpJdotRow;
		for (int jidx = 0; jidx < tmpJ.rows(); ++jidx) {
			if (tmpJ.row(jidx).norm() != 0) {
				tmpJRow = tmpJ.row(jidx),
					tmpJdotRow = tmpJdot.row(jidx);
				for (int pidx = 0; pidx < numParts; ++pidx) {
					tmpJ2(matLen, pidx) = tmpJRow[pidx];
					tmpJDot2(matLen, pidx) = tmpJdotRow[pidx];
				}//for each entry in row
				tmpCVal2[matLen] = tmpCVal[jidx];
				tmpCDotVal2[matLen] = tmpCDotVal[jidx];
				++matLen;
			}//if len != 0 means if non-zero entries in row
		}//for every row of old matrix
		if (matLen < numCnstrnts) {//if a zero-constraint row existed
			tmpJ.resize(matLen, numParts);
			tmpJdot.resize(matLen, numParts);
			tmpCVal.resize(matLen);
			tmpCDotVal.resize(matLen);

			tmpJ.block(0, 0, matLen, numParts) = tmpJ2.block(0, 0, matLen, numParts);
			tmpJdot.block(0, 0, matLen, numParts) = tmpJDot2.block(0, 0, matLen, numParts);

			tmpCVal = tmpCVal2.head(matLen);
			tmpCDotVal = tmpCDotVal2.head(matLen);

		}//if needing to rebuild martices to get rid of 0-contributing constraint rows
		 //now need to handle zero-particle cols

		 //set global vals
		J = (tmpJ);
		Jdot = (tmpJdot);
		q = (tmpq);
		qdot = (tmpqdot);
		Q = (tmpqdotdot);
		CVal = (tmpCVal);
		CDotVal = (tmpCDotVal);
	}//buildCnstrntStructJumper

    //build constraint structure for part 5 - constraint jumping - particle is only affected by contraint that is closest to where it is moving next timestep
	void mySystem::buildCnstrntStruct(bool multiCnstrnt){
		if (multiCnstrnt) { buildCnstrntStructJumper(); return; }
		//init J, Jdot matrices; lambda, C vecs, Cdot vecs
		int numParts = p.size()*3;
		int numCnstrnts = c.size();
		Eigen::MatrixXd tmpJ(numCnstrnts, numParts);			//1 col per particle's x,y,z (3n) , 1 row per constraint(m)
		Eigen::MatrixXd tmpJdot(numCnstrnts,numParts);			//1 col per particle's x,y,z (3n) , 1 row per constraint(m)
		tmpJ.setZero();
		tmpJdot.setZero();
		Eigen::VectorXd tmpCVal(numCnstrnts);							//c evaluated
		Eigen::VectorXd tmpCDotVal(numCnstrnts);						//c dot evaluated
		tmpCVal.setZero();
		tmpCDotVal.setZero();

		Eigen::VectorXd tmpq(numParts);							//q vector is position of all particles
		Eigen::VectorXd tmpqdot(numParts);						//qdot vector is velocity of all particles
		Eigen::VectorXd tmpqdotdot(numParts);							//Q vector is forces on each particle - qdoubledot will be WQ (matrix mult)
		tmpq.setZero();
		tmpqdot.setZero();
		tmpqdotdot.setZero();

		Eigen::Vector3d
			tmpVec4, tmpVec5, tmpVec42, tmpVec52;

		int nIdx = 0;							//index into part-related vals is nIdx, constraint-related vals is mIdx
		//first set up q, qdot, qdotdot vectors
		for (unsigned int idx = 0; idx < p.size(); ++idx) {
			//cout << "calc constraint for idx : " << idx << " part id : " << p[idx]->ID << " qIDX : " << (nIdx) << endl;
			tmpq.segment<3>(nIdx) = p[idx]->getPosition();
			tmpqdot.segment<3>(nIdx) = p[idx]->getVelocity();
			tmpqdotdot.segment<3>(nIdx) = p[idx]->getForceAcc();
			nIdx += 3;		//increment by 3 for each particle
		}
		int p1nIDX, p2nIDX;
		//next build constraint values for each constraint
		for (unsigned int cIdx = 0; cIdx < numCnstrnts; ++cIdx) {
			tmpVec4 = c[cIdx]->calcPartialCP1(p[c[cIdx]->p1Idx], p[c[cIdx]->p2Idx]);	//path constraint
			tmpVec5 = c[cIdx]->calcPartialCdotP1(p[c[cIdx]->p1Idx], p[c[cIdx]->p2Idx]);
			//tmpVec42 = c[cIdx]->calcPartialCP2(p[c[cIdx]->p1Idx], p[c[cIdx]->p2Idx]);	//path constraint
			//tmpVec52 = c[cIdx]->calcPartialCdotP2(p[c[cIdx]->p1Idx], p[c[cIdx]->p2Idx]);
			tmpVec42 = c[cIdx]->calcPartialCP1(p[c[cIdx]->p2Idx], p[c[cIdx]->p1Idx]);	//path constraint
			tmpVec52 = c[cIdx]->calcPartialCdotP1(p[c[cIdx]->p2Idx], p[c[cIdx]->p1Idx]);
			p1nIDX = 3 * c[cIdx]->p1Idx;
			p2nIDX = 3 * c[cIdx]->p2Idx;

			tmpCVal[cIdx] = c[cIdx]->ks * c[cIdx]->calcCVal(p[c[cIdx]->p1Idx], p[c[cIdx]->p2Idx]);
			tmpCDotVal[cIdx] = c[cIdx]->kd * c[cIdx]->calcCDotVal(p[c[cIdx]->p1Idx], p[c[cIdx]->p2Idx]);
			//SPD IMPLEMENTATION TODO
			//tmpCVal[cIdx] = c[cIdx].ks * c[cIdx].calcCValSPD(p[c[cIdx].p1Idx],p[c[cIdx].p2Idx], deltaT);		//.8 = ks
			//tmpCDotVal[cIdx] = c[cIdx].kd * c[cIdx].calcCDotValSPD(p[c[cIdx].p1Idx],p[c[cIdx].p2Idx], deltaT);//.8 = kd

			tmpJ.row(cIdx).segment<3>(p1nIDX) = tmpVec4;
			tmpJdot.row(cIdx).segment<3>(p1nIDX) = tmpVec5;
			tmpJ.row(cIdx).segment<3>(p2nIDX) = tmpVec42;
			tmpJdot.row(cIdx).segment<3>(p2nIDX) = tmpVec52;
		}

	
		//now need to remove from J matrix every row that is completely 0, otherwise JWJt will not be invertible, and need to remove equivalent row from Jdot, component from lambda vector and component from cval and cdotval vectors
		int matLen = 0;		//number of rows in new matrix
		Eigen::MatrixXd tmpJ2(numCnstrnts, numParts), tmpJDot2(numCnstrnts, numParts);			//will then shrink these once we have all 0-impact constraints removed
		Eigen::VectorXd tmpCVal2(numCnstrnts), tmpCDotVal2(numCnstrnts);
		Eigen::VectorXd tmpJRow, tmpJdotRow;
		for(int jidx = 0; jidx < tmpJ.rows(); ++jidx){
			if (tmpJ.row(jidx).norm() != 0){
				tmpJRow = tmpJ.row(jidx),
				tmpJdotRow = tmpJdot.row(jidx);
				for(int pidx = 0; pidx < numParts; ++pidx){
					tmpJ2(matLen,pidx) = tmpJRow[pidx];
					tmpJDot2(matLen,pidx) = tmpJdotRow[pidx];
				}//for each entry in row
				tmpCVal2[matLen] = tmpCVal[jidx];
				tmpCDotVal2[matLen] = tmpCDotVal[jidx];
				++matLen;
			}//if len != 0 means if non-zero entries in row
		}//for every row of old matrix
		if(matLen < numCnstrnts){//if a zero-constraint row existed
			tmpJ.resize(matLen, numParts);
			tmpJdot.resize(matLen, numParts);
			tmpCVal.resize(matLen);
			tmpCDotVal.resize(matLen);

			tmpJ.block(0, 0, matLen, numParts) = tmpJ2.block(0, 0, matLen, numParts);
			tmpJdot.block(0, 0, matLen, numParts) = tmpJDot2.block(0, 0, matLen, numParts);

			tmpCVal = tmpCVal2.head(matLen);
			tmpCDotVal = tmpCDotVal2.head(matLen);

		}//if needing to rebuild martices to get rid of 0-contributing constraint rows
		//now need to handle zero-particle cols

		//set global vals
		J = (tmpJ);
		Jdot = (tmpJdot);
		q = (tmpq);
		qdot = (tmpqdot);
		Q = (tmpqdotdot);
		CVal = (tmpCVal); 
		CDotVal = (tmpCDotVal);
		
	}//void buildConstraintStructure();

	void mySystem::handleTimeStep() {

	}


}//particleSystem namespace
