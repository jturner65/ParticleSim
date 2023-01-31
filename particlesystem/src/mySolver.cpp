#include "..\include\mySolver.h"
#include "..\include\myForce.h"

namespace particleSystem{
	unsigned int mySolver::ID_gen = 0;

	mySolver::mySolver(SolverType _t) :ID(++ID_gen), type(_t), lambda(2), invLam(0), lamHalf(0), hafMInvLam(0), oneMlamHalf(0)
	{ setLambda(lambda); setIntegratorType(_t); }
	mySolver::~mySolver() {}

	
	//USING single state and single statedot vectors
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
	Eigen::VectorXd mySolver::IntegrateGroundPerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot) {
		Eigen::VectorXd res(6), sDotDot(6); 
		sDotDot << _stateDot.segment<3>(3), 0, 0, 0;
		res = _state.segment<6>(0) + (deltaT *_stateDot.segment<6>(0)) + (.5 * deltaT * deltaT) * sDotDot;

		return res;
	}

	//explicit euler single time step evaluation
	Eigen::VectorXd mySolver::IntegrateExp_EPerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot) {
		Eigen::VectorXd res(6);
		res.setZero();
		res = _state.segment<6>(0) + (deltaT *_stateDot.segment<6>(0));
		return res;
	}

	//explicit euler single time step evaluation with force update - only called internally
	Eigen::VectorXd mySolver::IntegrateExp_EFrcUpd(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot) {
		Eigen::VectorXd res(9);
		res.setZero();
		res.segment<6>(0) = _state.segment<6>(0) + (deltaT *_stateDot.segment<6>(0));
		//add updated force here - recalc constraint frces TODO 
		res.segment<3>(6) = _stateDot.segment<3>(3);// +deltaT*(_stateDot.segment<3>(3) - _stateDot.segment<3>(9));		//finite diff for new force
		return res;
	}

	//midpoint method single time step evaulation - evaluate at deltaT/2, take the values found and use to re-evaluate
	Eigen::VectorXd mySolver::IntegrateMidpointPerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot) {
		Eigen::VectorXd tmpStateDot(6), res(6);
		Eigen::VectorXd deltaXhalf = IntegrateExp_EFrcUpd((deltaT *.5), _state, _stateDot);		//find new state at half-time
		
		//tmpStateDot = _stateDot;
		tmpStateDot.segment<6>(0) = deltaXhalf.segment<6>(3);				//copy over half-time vel
		res = IntegrateExp_EPerPart(deltaT, _state, tmpStateDot);	//x0 + h xdot1/2
		return res;
	}

	//state_t+1 = 2*state_t - state_t-1  + .5 * statedotdot_t * deltat*deltat
	//past state starts at segment loc 6+
	Eigen::VectorXd mySolver::IntegrateVerletPerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot) {
		//Eigen::VectorXd posNow = _state.segment<3>(0), posOld = _state.segment<3>(6),
		//	velNow = _state.segment<3>(3), velOld = _state.segment<3>(9),
		//	frcNow = _stateDot.segment<3>(3), frcOld = _stateDot.segment<3>(9);
		Eigen::VectorXd res(6);
		res.setZero();
		double delTsq = deltaT * deltaT;
		res.segment<3>(0) = _state.segment<3>(0) + ((_state.segment<3>(0) - _state.segment<3>(6)) * VERLET1mDAMP) + (delTsq * _stateDot.segment<3>(3));
		res.segment<3>(3) = _state.segment<3>(3) + ((_state.segment<3>(3) - _state.segment<3>(9)) * VERLET1mDAMP) + (delTsq * .5*(_stateDot.segment<3>(3) - _stateDot.segment<3>(9)));
	
		oldDelT = deltaT;//if delta t changes TODO
		return res;
	}

	Eigen::VectorXd mySolver::IntegrateRK3PerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot) {
		Eigen::VectorXd res(6); res.setZero();

		Eigen::VectorXd tmpVecK1(6), tmpVecK2(6), tmpVecK3(6); tmpVecK1.setZero(); tmpVecK2.setZero(); tmpVecK3.setZero();

		Eigen::VectorXd tmpVecState1 = IntegrateExp_EFrcUpd(deltaT, _state, _stateDot);
		tmpVecK1 << tmpVecState1.segment<6>(3);		//move resultant velocity and accel into xdot position

		Eigen::VectorXd tmpVecState2 = IntegrateExp_EFrcUpd((deltaT *.5), _state, tmpVecK1);
		tmpVecK2 << tmpVecState2.segment<6>(3);		//move resultant velocity and accel into xdot position

		Eigen::VectorXd tmpVecState3 = IntegrateExp_EFrcUpd(deltaT, _state, tmpVecK2);
		tmpVecK3 << tmpVecState3.segment<6>(3);		//move resultant velocity and accel into xdot position

		res = _state.segment<6>(0) + deltaT * ((tmpVecK1 + 4 * tmpVecK2 + tmpVecK3)/6.0);
		return res;
	}
	//RK4-solver
	Eigen::VectorXd mySolver::IntegrateRK4PerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot) {
		Eigen::VectorXd res(6); res.setZero();

		Eigen::VectorXd tmpVecK1(6), tmpVecK2(6), tmpVecK3(6), tmpVecK4(6);
		tmpVecK1.setZero(); tmpVecK2.setZero(); tmpVecK3.setZero();  tmpVecK4.setZero();

		Eigen::VectorXd tmpVecState1 = IntegrateExp_EFrcUpd(deltaT, _state, _stateDot);
		tmpVecK1 << tmpVecState1.segment<6>(3);		//move resultant velocity and accel into xdot position

		Eigen::VectorXd tmpVecState2 = IntegrateExp_EFrcUpd((deltaT *.5), _state, tmpVecK1);
		tmpVecK2 << tmpVecState2.segment<6>(3);		//move resultant velocity and accel into xdot position

		Eigen::VectorXd tmpVecState3 = IntegrateExp_EFrcUpd((deltaT *.5), _state, tmpVecK2);
		tmpVecK3 << tmpVecState3.segment<6>(3);		//move resultant velocity and accel into xdot position

		Eigen::VectorXd tmpVecState4 = IntegrateExp_EFrcUpd(deltaT, _state, tmpVecK3);
		tmpVecK4 << tmpVecState4.segment<6>(3);		//move resultant velocity and accel into xdot position
		
		//tmpVecK4.segment<3>(3) << tmpVecK3.segment<3>(3);			//move acceleration into vdot position

		res = _state.segment<6>(0) + deltaT * ((tmpVecK1 + 2 * (tmpVecK2 + tmpVecK3) + tmpVecK4) / 6.0);

		return res;
	}
	//general form as per Delin and Zhang
	Eigen::VectorXd mySolver::IntegrateRK4PerPartLambda(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot) {

		Eigen::VectorXd res(6);
		Eigen::VectorXd tmpVecK1(6), tmpVecK2(6), tmpVecK3(6), tmpVecK4(6), tmpVecK2a(6), tmpVecK3a(6);
		tmpVecK1.setZero(); tmpVecK2.setZero(); tmpVecK3.setZero(); tmpVecK4.setZero(); tmpVecK2a.setZero(); tmpVecK3a.setZero();

		Eigen::VectorXd tmpVecState1 = IntegrateExp_EFrcUpd(deltaT, _state, _stateDot);
		tmpVecK1 << tmpVecState1.segment<6>(3);		//move resultant velocity and accel into xdot position

		Eigen::VectorXd tmpVecState2 = IntegrateExp_EFrcUpd((deltaT *.5), _state, tmpVecK1);
		tmpVecK2 << tmpVecState2.segment<6>(3);		//move resultant velocity and accel into xdot position
		//mult by 2 because using hard 1/2 mult of timestep
		tmpVecK2a<< 2.0*((hafMInvLam * tmpVecK1) + (invLam * tmpVecK2));

		Eigen::VectorXd tmpVecState3 = IntegrateExp_EFrcUpd((deltaT *.5), _state, tmpVecK2a);
		tmpVecK3 << tmpVecState3.segment<6>(3);		//move resultant velocity and accel into xdot position

		tmpVecK3a.segment<6>(0) << ((oneMlamHalf * tmpVecK2.segment<6>(0)) + (lamHalf * tmpVecK3.segment<6>(0)));

		Eigen::VectorXd tmpVecState4 = IntegrateExp_EFrcUpd(deltaT, _state, tmpVecK3a);
		tmpVecK4 << tmpVecState4.segment<6>(3);		//move resultant velocity and accel into xdot position

		res = _state.segment<6>(0) + deltaT * ((tmpVecK1 + ((4 - lambda) * tmpVecK2) + (lambda * tmpVecK3) + tmpVecK4) / 6.0);

		return res;


	}
	//NOT WORKING - using conj grad solver in mySystem for implicit, but that only handles springs
	Eigen::VectorXd mySolver::IntegrateImp_EPerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot) {
		Eigen::VectorXd res(6);
		res.setZero();
		//Eigen::VectorXd tmpVec(6);
		//tmpVec.setZero();
		//Eigen::Vector3d delV(0, 0, 0);

		////Ynew = Y0 + delT * f(Ynew)
		////using taylor approx, f(Ynew) = f(Y0) + f'(Y0)*delY
		////delY = Ynew - Yold
		//// --> delY = (1/delT * I - f'(Y0))^-1 * f(Y0)
		//Eigen::Matrix3d I;  //3x3 identity
		//I.setIdentity();

		res = _state +(deltaT * _stateDot);
		//semi implicit - use t+1 vel
		res.segment<3>(0) = _state.segment<3>(0) + (deltaT * res.segment<3>(3));


													 //Eigen::VectorXd delY(2,Eigen::Vector3d(0,0,0));		// --> delY = (1/delT * I - f'(Y0))^-1 * f(Y0)
													 //delY[0] = 1.0/deltaT *
													 //delY[1] = 


													 //tmpVec[0] = _state[0] + (delY[0]);//pos				//xnew = x0 + delT * f(x0) * f'(xnew)
													 //tmpVec[1] = _state[1] + (delY[1]);					//vnew = v0 + delT * f(v0) * f'(vnew)

		return res;
	}

	Eigen::VectorXd mySolver::IntegrateTrapPerPart(double deltaT, const Eigen::Ref<const Eigen::VectorXd>& _state, const Eigen::Ref<const Eigen::VectorXd>& _stateDot) {
		Eigen::VectorXd res(6), resDot(9);
		res.setZero();
		//if(_perPart){}//{	std::cout<<"IntegrateTrapPerPart called per part\n"<<std::endl;	}
		//else	{		std::cout<<"IntegrateTrapPerPart called in loop\n"<<std::endl;	}
		resDot = IntegrateExp_EFrcUpd(deltaT, _state, _stateDot);//fwd euler w/frce upd
		//update force as well
		//trap - use 1/2 future and 1/2 past statedot
		res = _state.segment<6>(0) + (deltaT * .5 * (resDot.segment<6>(3) + _stateDot.segment<6>(0)));
		return res;
	}


}//particleSystem namespace
