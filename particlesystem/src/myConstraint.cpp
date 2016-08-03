#include "..\include\myConstraint.h"

namespace particleSystem{
	unsigned int myConstraint::ID_gen = 0;

	void myConstraint::setP1(std::shared_ptr<myParticle> _p, int _pIdx1){ 
        p1 = _p;
		p1ID = _p->ID;
		p1Idx = _pIdx1;
	}
	
	void myConstraint::setP2(std::shared_ptr<myParticle> _p, int _pIdx2){
        p2 = _p;
		p2ID = _p->ID;
		p2Idx = _pIdx2;
	}

}//namespace particleSystem
