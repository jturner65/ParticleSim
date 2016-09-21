#include "..\include\myConstraint.h"

namespace particleSystem{
	unsigned int myConstraint::ID_gen = 0;

	void myConstraint::setP1(std::shared_ptr<myParticle> _p, int _pIdx1){ 
        p1 = _p;
		p1ID = _p->ID;
		p1Idx = _pIdx1;
		if (p1ID != p1Idx) { std::cout << "p1Id : " << p1ID << " != p1Idx : " << p1Idx << "\n"; }
	}
	
	void myConstraint::setP2(std::shared_ptr<myParticle> _p, int _pIdx2){
        p2 = _p;
		p2ID = _p->ID;
		p2Idx = _pIdx2;
		if (p2ID != p2Idx) { std::cout << "p2ID : " << p2ID << " != p2Idx : " << p2Idx << "\n"; }
	}

}//namespace particleSystem
