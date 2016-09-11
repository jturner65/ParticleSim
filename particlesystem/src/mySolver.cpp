#include "..\include\mySolver.h"
#include "..\include\myForce.h"

namespace particleSystem{
	unsigned int mySolver::ID_gen = 0;
	//set function pointer based on passed type of solver
	void mySolver::setIntegratorType(SolverType type) {
		switch (type) {			
			case GROUND: {		integrator = &mySolver::IntegrateGroundPerPart; break; }
			case EXP_E: {		integrator = &mySolver::IntegrateExp_EPerPart; break; }
			case MIDPOINT: {		integrator = &mySolver::IntegrateMidpointPerPart; break; }
			case RK3: {		integrator = &mySolver::IntegrateRK3PerPart; break; }
			case RK4: {		integrator = &mySolver::IntegrateRK4PerPart; break; }
			case IMP_E: {		integrator = &mySolver::IntegrateImp_EPerPart; break; }
			case TRAP: {		integrator = &mySolver::IntegrateTrapPerPart; break; }
			case VERLET: {		integrator = &mySolver::IntegrateVerletPerPart; break; }
			case RK4_G: {       integrator = &mySolver::IntegrateRK4PerPartLambda; break; }
			default: {		integrator = &mySolver::IntegrateNone; break; }		//no solver specified
		}//switch
	}//setIntegratorType

	//derives the ground truth time step for a particle using pure physics and no external forces outside gravity
	vector<Eigen::Vector3d> mySolver::IntegrateGroundPerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot) {
		vector<Eigen::Vector3d> tmpVec(2, Eigen::Vector3d(0, 0, 0));
		tmpVec[0] = _state[0] + (_state[1] * deltaT) + (gravVec * (.5 * deltaT * deltaT));
		tmpVec[1] = _state[1] + (gravVec * deltaT);
		return tmpVec;
	}

	//explicit euler single time step evaluation
	vector<Eigen::Vector3d> mySolver::IntegrateExp_EPerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot) {
		vector<Eigen::Vector3d> tmpVec(2, Eigen::Vector3d(0, 0, 0));
		tmpVec[0] = _state[0] + (deltaT * _stateDot[0]);
		tmpVec[1] = _state[1] + (deltaT * _stateDot[1]);
		return tmpVec;
	}

	//midpoint method single time step evaulation - evaluate at deltaT/2, take the values found and use to re-evaluate
	vector<Eigen::Vector3d> mySolver::IntegrateMidpointPerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot) {
		vector<Eigen::Vector3d> tmpVec(2, Eigen::Vector3d(0, 0, 0));
		vector<Eigen::Vector3d> deltaXhalf = IntegrateExp_EPerPart((deltaT *.5), _state, _stateDot);
		vector<Eigen::Vector3d> tmpStateDot(2, Eigen::Vector3d(0, 0, 0));

		tmpStateDot[0] = deltaXhalf[1];			//new stateDot 0 term is v  @ t=.5 deltat, accel @ t = 0
		tmpStateDot[1] = _stateDot[1];			//deltaV is the same acceleration = _stateDot[1]
		tmpVec = IntegrateExp_EPerPart(deltaT, _state, tmpStateDot);	//x0 + h xdot1/2
		return tmpVec;
	}

	//Verlet method single time step evaulation (using old values to approximate)
	//_state[0] : cur pos
	//_state[1] : cur vel
	//_state[2] : last pos
	//_state[3] : last vel
	//_stateDot[0] : cur vel
	//_stateDot[1] : cur f
	//_stateDot[2] : last vel
	//_stateDot[3] : last f

	//state_t+1 = 2*state_t - state_t-1  + .5 * statedotdot_t * deltat*deltat
	vector<Eigen::Vector3d> mySolver::IntegrateVerletPerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot) {
		//vector<Eigen::Vector3d> tmpMPRes = IntegrateExp_EPerPart((deltaT *.5), _state, _stateDot);     //for midpoint velocity - tmpMPRes[1]
		vector<Eigen::Vector3d> tmpVec(2, Eigen::Vector3d(0, 0, 0));

		//tmpVec[0] = _state[0] + deltaT * tmpMPRes[1];        //velocity @ .5 timestep

		//verlet : pos_t1 = 2*pos_t0 - pos_t-1 + delt^2 * accel
		tmpVec[0] = _state[0] + ((_state[0] - _state[2]) * VERLET1mDAMP) + (deltaT * deltaT * _stateDot[1]);           //verlet without velocity    
		tmpVec[1] = _state[1] + ((_state[1] - _state[3]) * VERLET1mDAMP) + (deltaT * deltaT *(.5* (_stateDot[1] - _stateDot[3])));           //verlet without velocity    

		//tmpVec[1] = (tmpVec[0] - _state[0])/(deltaT) + (deltaT * _stateDot[1]);         //need velocity for collisions - mean value thm + fwd euler acc? not very accurate


		oldDelT = deltaT;//if delta t changes TODO
		return tmpVec;
	}

	vector<Eigen::Vector3d> mySolver::IntegrateRK3PerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot) {
		//if(_perPart){}//{	std::cout<<"IntegrateRK3PerPart called per part\n"<<std::endl;	}
		//else	{		std::cout<<"IntegrateRK3PerPart called in loop\n"<<std::endl;	}
		vector<Eigen::Vector3d> tmpVec(2, Eigen::Vector3d(0, 0, 0));
		vector<Eigen::Vector3d> tmpVecK1(2, Eigen::Vector3d(0, 0, 0));
		vector<Eigen::Vector3d> tmpVecK2(2, Eigen::Vector3d(0, 0, 0));
		vector<Eigen::Vector3d> tmpVecK3(2, Eigen::Vector3d(0, 0, 0));
		//vector<Eigen::Vector3d> tmpVecState1(2, Eigen::Vector3d(0,0,0));
		//vector<Eigen::Vector3d> tmpVecState2(2, Eigen::Vector3d(0,0,0));
		//vector<Eigen::Vector3d> tmpVecState3(2, Eigen::Vector3d(0,0,0));
		/**
		yi+1 = yi + 1/6 ( k1 + 4 k2 + k3 ),:k1 = h f(xi, yi),k2 = h f(xi + h / 2, yi + k1 / 2 ),k3 = h f(xi + h, yi - k1 + 2 k2 ),
		*/

		//vector<Eigen::Vector3d> tmpVecState1 = IntegrateExp_EPerPart(deltaT, _state, _stateDot);
		tmpVecK1[0] = _state[1];		//move resultant velocity into xdot position
		tmpVecK1[1] = _stateDot[1];			//move acceleration into vdot position

		vector<Eigen::Vector3d> tmpVecState2 = IntegrateExp_EPerPart((deltaT *.5), _state, tmpVecK1);
		tmpVecK2[0] = tmpVecState2[1];		//move resultant velocity into xdot position
		tmpVecK2[1] = tmpVecK1[1];			//move acceleration into vdot position

		vector<Eigen::Vector3d> tmpVecState3 = IntegrateExp_EPerPart(deltaT, _state, tmpVecK2);
		tmpVecK3[0] = tmpVecState3[1];
		tmpVecK3[1] = tmpVecK2[1];			//tmpVecK3 should just be delta part of exp euler evaluation

		tmpVec[0] = _state[0] + deltaT * ((tmpVecK1[0] + 4 * tmpVecK2[0] + tmpVecK3[0]) / 6.0);
		tmpVec[1] = _state[1] + deltaT * ((tmpVecK1[1] + 4 * tmpVecK2[1] + tmpVecK3[1]) / 6.0);

		return tmpVec;
	}
	//RK4-solver
	vector<Eigen::Vector3d> mySolver::IntegrateRK4PerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot) {
		vector<Eigen::Vector3d> tmpVec(2, Eigen::Vector3d(0, 0, 0));
		vector<Eigen::Vector3d> tmpVecK1(2, Eigen::Vector3d(0, 0, 0));
		vector<Eigen::Vector3d> tmpVecK2(2, Eigen::Vector3d(0, 0, 0));
		vector<Eigen::Vector3d> tmpVecK3(2, Eigen::Vector3d(0, 0, 0));
		vector<Eigen::Vector3d> tmpVecK4(2, Eigen::Vector3d(0, 0, 0));

		//vector<Eigen::Vector3d> tmpVecState1 = IntegrateExp_EPerPart(deltaT, _state, _stateDot);
		tmpVecK1[0] = _state[1];		//start state 
		tmpVecK1[1] = _stateDot[1];			//move acceleration into vdot position

		vector<Eigen::Vector3d> tmpVecState2 = IntegrateExp_EPerPart((deltaT *.5), _state, tmpVecK1);
		tmpVecK2[0] = tmpVecState2[1];		            //move resultant velocity into xdot position
		tmpVecK2[1] = tmpVecK1[1];			//move acceleration into vdot position

		vector<Eigen::Vector3d> tmpVecState3 = IntegrateExp_EPerPart((deltaT *.5), _state, tmpVecK2);
		tmpVecK3[0] = tmpVecState3[1];
		tmpVecK3[1] = tmpVecK2[1];			//tmpVecK3 should just be delta part of exp euler evaluation

		vector<Eigen::Vector3d> tmpVecState4 = IntegrateExp_EPerPart(deltaT, _state, tmpVecK3);
		tmpVecK4[0] = tmpVecState4[1];
		tmpVecK4[1] = tmpVecK3[1];			//tmpVecK3 should just be delta part of exp euler evaluation

		tmpVec[0] = _state[0] + deltaT * ((tmpVecK1[0] + 2 * (tmpVecK2[0] + tmpVecK3[0]) + tmpVecK4[0]) / 6.0);
		tmpVec[1] = _state[1] + deltaT * ((tmpVecK1[1] + 2 * (tmpVecK2[1] + tmpVecK3[1]) + tmpVecK4[1]) / 6.0);

		return tmpVec;
	}
	//general form as per Delin and Zhang
	vector<Eigen::Vector3d> mySolver::IntegrateRK4PerPartLambda(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot) {
		vector<Eigen::Vector3d> tmpVec(2, Eigen::Vector3d(0, 0, 0));
		vector<Eigen::Vector3d> tmpVecK1(2, Eigen::Vector3d(0, 0, 0));
		vector<Eigen::Vector3d> tmpVecK2(2, Eigen::Vector3d(0, 0, 0));
		vector<Eigen::Vector3d> tmpVecK3(2, Eigen::Vector3d(0, 0, 0));
		vector<Eigen::Vector3d> tmpVecK2a(2, Eigen::Vector3d(0, 0, 0));
		vector<Eigen::Vector3d> tmpVecK3a(2, Eigen::Vector3d(0, 0, 0));
		vector<Eigen::Vector3d> tmpVecK4(2, Eigen::Vector3d(0, 0, 0));

		tmpVecK1[0] = _state[1];		//move resultant velocity into xdot position
		tmpVecK1[1] = _stateDot[1];			//move acceleration into vdot position

		vector<Eigen::Vector3d> tmpVecState2 = IntegrateExp_EPerPart((deltaT *.5), _state, tmpVecK1);
		tmpVecK2[0] = tmpVecState2[1];		//move resultant velocity into xdot position - general form uses 1 and 2
		tmpVecK2[1] = tmpVecK1[1];			//move acceleration into vdot position

		tmpVecK2a[0] = (((.5 - (1.0 / lambda)) * _state[1]) + ((1.0 / lambda) * tmpVecState2[1]));		//move resultant velocity into xdot position - general form uses 1 and 2
		tmpVecK2a[1] = tmpVecK1[1];			//move acceleration into vdot position

		vector<Eigen::Vector3d> tmpVecState3 = IntegrateExp_EPerPart((deltaT *.5), _state, tmpVecK2a);
		tmpVecK3[0] = tmpVecState3[1];
		tmpVecK3[1] = tmpVecK2[1];			//tmpVecK3 should just be delta part of exp euler evaluation

		tmpVecK3a[0] = (((1 - (lambda / 2.0)) * tmpVecState2[1]) + ((lambda / 2.0) * tmpVecState3[1]));
		tmpVecK3a[1] = tmpVecK2a[1];			//tmpVecK3 should just be delta part of exp euler evaluation

		vector<Eigen::Vector3d> tmpVecState4 = IntegrateExp_EPerPart(deltaT, _state, tmpVecK3a);
		tmpVecK4[0] = tmpVecState4[1];
		tmpVecK4[1] = tmpVecK3[1];			//tmpVecK3 should just be delta part of exp euler evaluation

		tmpVec[0] = _state[0] + deltaT * ((tmpVecK1[0] + ((4 - lambda) * tmpVecK2[0]) + (lambda * tmpVecK3[0]) + tmpVecK4[0]) / 6.0);
		tmpVec[1] = _state[1] + deltaT * ((tmpVecK1[1] + ((4 - lambda) * tmpVecK2[1]) + (lambda * tmpVecK3[1]) + tmpVecK4[1]) / 6.0);

		return tmpVec;
	}
	//NOT WORKING - using conj grad solver in mySystem for implicit, but this only handles springs
	vector<Eigen::Vector3d> mySolver::IntegrateImp_EPerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot) {
		vector<Eigen::Vector3d> tmpVec(2, Eigen::Vector3d(0, 0, 0));
		Eigen::Vector3d delV(0, 0, 0);

		//Ynew = Y0 + delT * f(Ynew)
		//using taylor approx, f(Ynew) = f(Y0) + f'(Y0)*delY
		//delY = Ynew - Yold
		// --> delY = (1/delT * I - f'(Y0))^-1 * f(Y0)
		Eigen::Matrix3d I;  //3x3 identity
		I.setIdentity();

		//velocity 
		tmpVec[1] = _state[1] + (deltaT * _stateDot[1]);// + (deltaT * );//v : _stateDot[1] is f(v0)  we want f(v1) = f(v0) + delV * f'(v0) == delV = (1/delT * I - f'(v0))^-1 * f(v0)

		//have Vnew to calc new position
		tmpVec[0] = _state[0] + (deltaT * tmpVec[1]);//pos			//tmpVec[1] = v(t+dt)


		//vector<Eigen::Vector3d> delY(2,Eigen::Vector3d(0,0,0));		// --> delY = (1/delT * I - f'(Y0))^-1 * f(Y0)
		//delY[0] = 1.0/deltaT *
		//delY[1] = 


		//tmpVec[0] = _state[0] + (delY[0]);//pos				//xnew = x0 + delT * f(x0) * f'(xnew)
		//tmpVec[1] = _state[1] + (delY[1]);					//vnew = v0 + delT * f(v0) * f'(vnew)

		return tmpVec;
	}
	vector<Eigen::Vector3d> mySolver::IntegrateTrapPerPart(double deltaT, vector<Eigen::Vector3d>& _state, vector<Eigen::Vector3d>& _stateDot) {
		vector<Eigen::Vector3d> tmpVec(2, Eigen::Vector3d(0, 0, 0));
		//if(_perPart){}//{	std::cout<<"IntegrateTrapPerPart called per part\n"<<std::endl;	}
		//else	{		std::cout<<"IntegrateTrapPerPart called in loop\n"<<std::endl;	}

		tmpVec[1] = _state[1] + (deltaT * _stateDot[1]);		//assuming const accelerations allow use of _statDot[1] - otherwise need to calculate for f(v(t+dt))
		tmpVec[0] = _state[0] + (deltaT * ((.5*tmpVec[1]) + (.5 * _state[1])));			//tmpVec[1] = v(t+dt)

		return tmpVec;
	}
	
}//particleSystem namespace
