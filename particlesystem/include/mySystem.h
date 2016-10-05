#ifndef mySystem_h
#define mySystem_h

#include <vector>

#include <Eigen/Dense>

#include "myParticle.h"
#include "mySolver.h"
#include "myInvPend.h"
#include "myForce.h"
#include "myConstraint.h"
#include "myCollider.h"
#include "myFluidBox.h"
#include "mySpring.h"

using namespace std;

namespace particleSystem{
	//should include structure to hold particles, forces, time
	class mySystem {
	public:

		mySystem(string _name, double _delT, int numP, int numF, int numC);
		~mySystem();

		inline void buildFluidBox(int numCellX, int numCellY, int numCellZ, double _diffusion, double _viscosity, const Eigen::Ref<const Eigen::Vector3d>& _ctrLoc, const Eigen::Ref<const Eigen::Vector3d>& cellSize) {
			fluidBox = make_shared<myFluidBox>(numCellX, numCellY, numCellZ, _diffusion, _viscosity, deltaT, cellSize);
			fluidBox->setCenter (Eigen::Vector3d(_ctrLoc));
		}//buildFluidBox
		//initialize system variables - run after any particles, constraints, colliders have been added or changed
		void initMassSystem();

		//overridden by individual systems based on configuration
		virtual void handlePause() {};
		virtual void handleTimeStep();

		inline void addParticle(double _mass, const Eigen::Ref<const Eigen::Vector3d>&  _pos, SolverType _solv) {
			p.push_back(std::make_shared<myParticle>(1.0, _pos, Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), _solv));
		}

		inline void addParticle(double _mass, const Eigen::Ref<const Eigen::Vector3d>&  _pos, const Eigen::Ref<const Eigen::Vector3d>&  _vel, SolverType _solv) {
			p.push_back(std::make_shared<myParticle>(1.0, _pos, _vel, Eigen::Vector3d(0, 0, 0), _solv));
		}

		void applyForcesToSystem();
		void applySpringForcesToSystem();

		void addForcesToTinkerToys(double fMult, bool msDragged, const Eigen::Ref<const Eigen::Vector3d>& msDragDif);
		void addForcesToMassSpringCtrl(double fMult, bool msDragged, const Eigen::Ref<const Eigen::Vector3d>& msDragDif);
		void addDraggingForceTT(Eigen::Vector3d& force, double mult);
		void addShakeForceToFluid(double fMult, const Eigen::Ref<const Eigen::Vector3d>& msdrgVal0, const Eigen::Ref<const Eigen::Vector3d>& msdrgVal1);
		void addControlMS2(double& tsCntr, int numMSMsclSpr);

		void invokeSolverDerivEval();

		inline void setNumParts(int _rs) { p.resize(_rs); }
		inline void setNumForces(int _rs) { f.resize(_rs); }
		inline void setNumConstraints(int _rs) { c.resize(_rs); }

		inline void setDeltaT(double _delT) { deltaT = _delT; }

		inline string buildDblStr(double val, const char* fmt) {                                          //build label to include current slider value
			char buf[256];
			stringstream ss;
			sprintf(buf, fmt, val);
			ss << buf;
			return ss.str();
		}

		//builds paren surrounded string of eigen vec3d vals
		inline std::string evec3dToStr(const Eigen::Ref<const Eigen::Vector3d>& in, const char* fmt = "%.4f") {
			stringstream ss;
			ss << "(" << buildDblStr(in(0), fmt) << "," << buildDblStr(in(1), fmt) << "," << buildDblStr(in(2), fmt) << ")";
			return ss.str();
		}

		inline double dist2Parts(int p1idx, int p2idx) { return (p[p1idx]->getPosition() - p[p2idx]->getPosition()).norm(); }

		//void buildRhTrHdrn(Eigen::Vector3d& sLoc); //rhombic triacontahedron (30 sided die) - for mass spring motion project
		//void buildMsSprMtn2(Eigen::Vector3d& sLoc); //some other mass-spring object
		//void buildFlatCollider(double krest, double muFrict, double x, double y, double z, const Eigen::Ref<const Eigen::Vector3d>&std::string& name, const Eigen::Ref<const Eigen::Vector3d>& drLoc, const Eigen::Ref<const Eigen::Vector3d>& gndLoc);
		void buildGndCollider(double krest, double muFrict, double x, double y, double z, std::string& name, const Eigen::Ref<const Eigen::Vector3d>& drLoc);// , const Eigen::Ref<const Eigen::Vector3d>& gndLoc);
		void buildGlobeCollider(double krest, double muFrict, double rad, double distFromGlb);

		void mySystem::buildDefForces() {		_buildBaseDefForces(0);	}//no drag
		void mySystem::buildDefTTForces() {		_buildBaseDefForces(dragTTFrcCoef);	}
		void mySystem::buildDefMSForces() {		_buildBaseDefForces(dragMSFrcCoef);	}
		void mySystem::buildDefFluidForces() {	_buildBaseDefForces(dragFluidFrcCoef);	}

	private : 
		void _buildBaseDefForces(double dCoeff);
		void buildCnstrntStructJumper();

		//inline void calcFeedbackTerm(unsigned int cIdx) {
		//	feedBack(cIdx) = c[cIdx]->ks * c[cIdx]->calcCVal(p[c[cIdx]->p1Idx], p[c[cIdx]->p2Idx]) + c[cIdx]->kd * c[cIdx]->calcCDotVal(p[c[cIdx]->p1Idx], p[c[cIdx]->p2Idx]);

		//}
		inline void calcFBTermSPD() {
			//Eigen::MatrixXd  invM = (mSkel->getMassMatrix() + mKd * mTimestep).inverse();
			//Eigen::VectorXd p = -mKp * (_dof + _dofVel * mTimestep - mDesiredDofs);                    //mDesiredDofs is where we want to end up
			//Eigen::VectorXd d = -mKd * _dofVel;
			//Eigen::VectorXd qddot = invM * (-mSkel->getCombinedVector() + p + d + mConstrForces);
			//mTorques = p + d - mKd * qddot * mTimestep;

			//Eigen::VectorXd q = mBiped->getPositions();
			//Eigen::VectorXd dq = mBiped->getVelocities();

			//Eigen::MatrixXd invM = (mBiped->getMassMatrix() + mKd * mBiped->getTimeStep()).inverse();
			//Eigen::VectorXd p = -mKp * (q + dq * mBiped->getTimeStep() - mTargetPositions);
			//Eigen::VectorXd d = -mKd * dq;
			//Eigen::VectorXd qddot = invM * (-mBiped->getCoriolisAndGravityForces() + p + d + mBiped->getConstraintForces());

			//feedBack = p + d - mKd * qddot * mBiped->getTimeStep();
		}

	public :

		void buildCnstrntStruct(bool multiCnstrnt);

		void buildInvPendChain(Eigen::Vector3d& sLoc, int numParts);
		void buildRollerCoasterConstraints(int id, double rad);
		void buildAndSetCnstrnts(int part1IDX, int part2IDX, double rad, Eigen::Vector3d& center);

		//void buildAndSetCnstrntVals(int partIDX, double rad, Eigen::Vector3d& center);
		//void buildAndSetNonPathCnstrntVals(int part1IDX, int part2IDX, Eigen::Vector3d& center);

		void resetCtrPointMassSpr();
		inline void reInitAllParts() { for (unsigned int i = 0; i<p.size(); ++i) { p[i]->reInitPartComps(); } }//resetCtrPointMassSpr();}

		inline int getParticleIDXByID(int id) { for (unsigned int i = 0; i<p.size(); ++i) { if (p[i]->ID == id) return i; } return -1; }
		inline int getCnstIDXByID(int id) { for (unsigned int i = 0; i<c.size(); ++i) { if (c[i]->ID == id) return i; } return -1; }

		Eigen::Vector3d calcAndSetCOM(unsigned int stIdx, unsigned int endIdx);               //find com of particles in sim in container, from idx stIdx to idx endIdx
		Eigen::Vector3d calcDragForce(Eigen::Vector3d& sValP);

		void calcAnkleForce();
		double calcAndApplyAnkleForce(int cidx);                //calculate and apply appropriate ankle forces to counteract forces on constrained particle, return kp (with kd = 1/5 * deltat * kp)
		void calcFluidForceForAllParticles();
		void dispCenterPart(Eigen::Vector3d& dir);

		//satisfy constraints positionally
		void satisfyPosConstraints();
		//solve implicit euler for mass springs here, instead of in solver
		void solveImpEuler4MassSpring();
		//conjugate gradient solver
		Eigen::VectorXd calcConjGrad(const Eigen::Ref<const Eigen::VectorXd>& b, const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::VectorXd>& f0);

		bool handleSnowGlobeTimeStep(double fmult, bool msDragged, const Eigen::Ref<const Eigen::Vector3d>& msDragVal0, const Eigen::Ref<const Eigen::Vector3d>& msDragVal1) {
			bool res = msDragged;
			//apply forces
			applyForcesToSystem();
			//handle fluid force stuff - res might change, is passed by ref
			if ((res) || ((shakeVal).norm() > 0.001)) {
				addShakeForceToFluid(fmult, msDragVal0, msDragVal1); 		res = false;			//setting false keeps force from compounding
			}
			fluidBox->myFluidSphereTimeStep();					//timestep snowglobe	
																//fluidBox->resetOldVals();
			calcFluidForceForAllParticles();				//calculate fluid forces
			shakeVal = Eigen::Vector3d(0, 0, 0);
			handlePartCldrCollision();									//checks for, and handles, collisions for all particles against colliders
			//handlePartPartRepelCollision();								//handle particle-particle interraction for snowglobe by applying repelling forces - assume particles are too small to collide
			invokeSolverDerivEval();
			return res;
		}

		void handleMassSpringTimeStep(bool addCntrl, bool msClicked, bool msDragged, int numMSMsclSpr, const Eigen::Ref<const Eigen::Vector3d>& msDragVal) {
			//reinitialize all particles color and mass
			reInitAllParts();
			if (addCntrl) {//handle animation-driven control - change addCntrl to flags var
				if (flags[mySystem::mass2HasHitGrnd]) {	addControlMS2(tsCntr, numMSMsclSpr);}				//when hit ground add motion control - derive appropriate spring lengths for desired motion				
			}
			else {
				if (msClicked) { shakeVal.setZero(); }
				if (msDragged) { addForcesToMassSpringCtrl(5, msDragged, msDragVal);	}
			}
			applyForcesToSystem();
			applySpringForcesToSystem();
			bool colTest = handlePartCldrCollision();			//checks for, and handles, collisions for all particles to see if they hit the ground
			if (colTest && flags[mass2ChkHasHitGrnd]) {			//check if ought to hit ground - only check 1 time, latches off check
				flags[mass2HasHitGrnd] = true;
				flags[mass2ChkHasHitGrnd] = false;
			}			
			satisfyPosConstraints();							//solve for positional constraints
			solveImpEuler4MassSpring();							//invoke implicit solver using cojugate gradient method
		}//handleMassSpringTimeStep

		void handleTTTimeStepRepel(double fmult, bool msDragged, bool jumpCnstrnt, const Eigen::Ref<const Eigen::Vector3d>& msDragVal) {
			applyForcesToSystem();
			if (msDragged) {		addForcesToTinkerToys(fmult, msDragged, msDragVal);		}	//adding dragged forces to system here
			handlePartPartRepelCollision();														//handle particle-particle interraction, with repeling force if close
			buildCnstrntStruct(jumpCnstrnt);													//handle constraints here - rebuild constraint structure
			calcConstraintForces();																//handle constraints here - calculate constraint forces and matrices
			applyConstraintForcesToSystem();													//handle collisions from contraint enforcement
			invokeSolverDerivEval();															//invoke derivative handler
		}//handleTTTimeStepRepel

		void handleTTTimeStep(double fmult, bool msDragged, bool jumpCnstrnt, const Eigen::Ref<const Eigen::Vector3d>& msDragVal) {
			applyForcesToSystem();
			if (msDragged) { addForcesToTinkerToys(fmult, msDragged, msDragVal); }				//adding dragged forces to system here
			handlePartPartCollision();															//handle particle-particle interraction, with repeling force if close
			buildCnstrntStruct(jumpCnstrnt);													//handle constraints here - rebuild constraint structure
			calcConstraintForces();																//handle constraints here - calculate constraint forces and matrices
			applyConstraintForcesToSystem();													//handle collisions from contraint enforcement
			invokeSolverDerivEval();															//invoke derivative handler
		}//handleTTTimeStep

		bool handleInvPendTimeStep(double fmult, bool hasFluidGlobe, bool jumpCnstrnt, bool msDragged, bool calcCOM, const Eigen::Ref<const Eigen::Vector3d>& msDragVal0, const Eigen::Ref<const Eigen::Vector3d>& msDragVal1) {
			bool res = msDragged;
			//constrained tinker toy code here
			applyForcesToSystem();
			if (hasFluidGlobe) {	//inv pend within fluid
				if ((res) || ((shakeVal).norm() > 0.001)) {
					addShakeForceToFluid(fmult, msDragVal0, msDragVal1); 		res = false;			//setting false keeps force from compounding
				}
				fluidBox->myFluidBoxTimeStep();					//timestep snowglobe	

				calcFluidForceForAllParticles(); 
				
				//handleGlobeTimeStep(fmult, res, msDragVal0, msDragVal1);
				handlePartCldrCollision();//checks for, and handles, collisions for all particles to see if they hit the ground
				handlePartPartCollision();
			}
			else {//#INV_PEND    - no fluid just add force        
				if (msDragged) {
					Eigen::Vector3d msDiff = msDragVal0 - msDragVal1;
					addForcesToTinkerToys(fmult, msDragged, msDiff);
				}
			}
			//calcAnkleForce();													//calc ankle forces based on derived forces

			//TODO make all ankle forces and constraint structures calculated per seaweed strand (make a class)
			for (unsigned int i = 0; i < invPend.size(); ++i) {
				invPend[i]->calcCnstrntFrc();
				invPend[i]->applyConstraintForcesToSystem();
			}
			if (calcCOM) {
				for (unsigned int i = 0; i < invPend.size(); ++i) {
					invPend[i]->calcAndSetCOM();
				}
			}

			invokeSolverDerivEval();											//invoke derivative handler

			return res;
		}

		void handlePartPartCollision();
		void handlePartPartRepelCollision();
		bool checkParticlePairForCollision(std::shared_ptr<myParticle> partA, std::shared_ptr<myParticle> partB, double dist);

	
		bool handlePartCldrCollision();
		//void handlePlanarCollision(myCollider& col, std::shared_ptr<myParticle> part, int colType);
		//void handleSphereCollision(myCollider& col, std::shared_ptr<myParticle> part, int colType);

		bool handlePauseDrawClick(int& clickOnPartIDX, const Eigen::Ref<const Eigen::Vector3d>& mClickLoc, bool idxIs5);
		bool handlePauseDrawRel(int& clickOnPartIDX, const Eigen::Ref<const Eigen::Vector3d>& mClickRel, bool idxIs5);

		void applyForceVecToAllPartsMsSpr(const Eigen::Ref<const Eigen::Vector3d>& frc);

		//apply all resultant constraint forces to system (qhat)
		void applyConstraintForcesToSystem() {
			unsigned int pSize = p.size();
			for (unsigned int pidx = 0; pidx < pSize; ++pidx) {
				p[pidx]->applyForce(Eigen::Vector3d(Qhat[pidx * 3], Qhat[(pidx * 3) + 1], Qhat[(pidx * 3) + 2]));
			}
		}//applyConstraintForcesToSystem

		void calcMassAndWMat();

		void calcConstraintForces();

		friend ostream& operator<<(ostream& out, const mySystem& sys) {
			out << "Delta t: " << sys.deltaT;// << " use particle-specific solvers? " << sys.usePartSolver << " system solver type : " << SolverType2str[sys.solveType];
			out << " Size of particle vector : " << sys.p.size();
			if (sys.p.size() >0) { out << endl; for (unsigned int idx = 0; idx < sys.p.size(); ++idx) { out << "\t" << sys.p[idx]; } }//if it has particles			
			out << " Size of force vector : " << sys.f.size();
			if (sys.f.size() >0) { out << endl; for (unsigned int idx = 0; idx < sys.f.size(); ++idx) { out << "\t" << sys.f[idx]; } }//if it has forces						
			out << " Size of collider vector : " << sys.colliders.size() << endl;
			if (sys.colliders.size() >0) { out << endl; for (unsigned int idx = 0; idx < sys.colliders.size(); ++idx) { out << "\t" << sys.colliders[idx]; } }//if it has constraints			
			out << " Size of constraint vector : " << sys.c.size() << endl;
			if (sys.c.size() >0) { out << endl; for (unsigned int idx = 0; idx < sys.c.size(); ++idx) { out << "\t" << sys.c[idx]; } }//if it has constraints			
			return out;
		}//op<<
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	public:
		static unsigned int ID_gen;

		int ID;
		string name;

		//Eigen::Vector3d ground;
		Eigen::Vector3d shakeVal;
		Eigen::Vector3d partCOM;

		vector<double> kpAra;

		vector<bool> flags;
		static const int debug = 0;								//init to false
		static const int pauseSim = 1;							//init to false
		static const int mass2HasHitGrnd = 2;					//init to true, eventually control via UI
		static const int mass2ChkHasHitGrnd = 3;				//init to true, eventually control via UI
		static const int buildDragCnstrnt = 4;					//init to false
		static const int canBuildDragPathCnstrnt = 5;			//init to false
		static const int showVel = 6;							//init to true; eventually control via UI
		static const int useMassMat = 7;						//whether we need to use the mass/w matrix (constraint solver) - init to true, turn off for snow globe which is huge

		static const int numFlags = 8;

		double deltaT, tsCntr;							//timestep
		std::shared_ptr<mySolver> solver;

		vector<std::shared_ptr<myParticle>> p;					//array holding particles
		vector<std::shared_ptr<myForce>> f;						//array holding all active forces in system
		vector<std::shared_ptr<myConstraint>> c;				//array holding all active constraints on system
		vector<std::shared_ptr<myCollider>> colliders;			//array holding the collision types the system will support
		vector<std::shared_ptr<myInvPend>> invPend;				//inverted pendulum ara/seaaweed ara
		vector<std::shared_ptr<mySpring>> spr;

		//constrained bead on wire matrices and vectors	n = # parts, m = num constraints
		Eigen::MatrixXd J,								//1 col per particle's x,y,z (3n) , 1 row per contraint(m)
			Jdot,							//1 col per particle's x,y,z (3n) , 1 row per contraint(m)
			W,								//diagonal inverse mass matrix, 1 3x3 diag entry per particle, rest 0 - only stores non-zero values in rows of sparse vectors (3n x 3n)
		//	spdInvM,						//inverse of mass matrix + kd*delT
			M;								//diagonal mass matrix, 1 3x3 diag entry per particle, rest 0 - only stores non-zero values in rows of sparse vectors (3n x 3n)
		Eigen::VectorXd q,										//q vector is position of all particles
			qdot,									//velocity of all particles
			Q,										//Q is vector of forces on each particle, from force acc
			feedBack,								//feedback term - pd control/spd control
			Qhat;// ,									//Qhat vector is constrain forces to apply to each particle - this is what is -applied- (result of constraint calculations) to force acc
		//	lambda;								//lambda is c.size() # of lagrangian multipliers, 1 per constraint

		std::shared_ptr<myFluidBox> fluidBox;								//fluid box to simulate snowglobe accelerations

	};//class mySystem

}//partsys
#endif