#ifndef myForce_h
#define myForce_h

#include <vector>
#include <string>
#include <iostream>
#include <Eigen/Dense>
#include "myParticle.h"
#include <memory>

using namespace std;

namespace particleSystem{

	class myForce {
	public:
		myForce(void) :ID(++ID_gen), name(""), constVal(), constVal2(), constVec(0, 0, 0), ftype(S_SCALAR) {}
		myForce(string _n, const Eigen::Vector3d& _G, ForceType _t) :ID(++ID_gen), name(_n), constVal(0), constVal2(0), constVec(_G), ftype(_t) {}		//
		myForce(string _n, double _k, ForceType _t) :ID(++ID_gen), name(_n), constVal(_k), constVal2(0), constVec(0, 0, 0), ftype(_t) {}		//if drag, needs to be negative constant value
		myForce(string _n, double _k) :ID(-1), name(_n), constVal(_k * (_k>0) ? 1 : -1), constVal2(0), constVec(0, 0, 0), ftype((_k>0) ? REPL : ATTR) {}						//repulsive force - dummy force, no incrementing id 
		myForce(string _n, double _k1, double _k2) :ID(++ID_gen), name(_n), constVal(_k1), constVal2(_k2), constVec(0, 0, 0), ftype(DAMPSPRING) {}		//2 constants are always SPRING
		myForce(string _n, double _k1, double _k2, ForceType _t) :ID(++ID_gen), name(_n), constVal(_k1), constVal2(_k2), constVec(0, 0, 0), ftype(_t) {}		//torque-result force
		myForce(const myForce& other) :ID(other.ID), name(other.name), constVal(other.constVal), constVal2(other.constVal2), constVec(other.constVec), ftype(other.ftype) {}
		~myForce(void) {}

		//pos 0 in vector is point 0, pos 1 is point 1
		static vector<Eigen::Vector3d> calcForceOnParticle(std::shared_ptr<myParticle> _p1, std::shared_ptr<myParticle> _p2, double d, std::shared_ptr<myForce> force);//sets quantities in particle relevant to this force application - adds force result to forceAcc

		//ForceType getFtype(){return ftype;}

		friend ostream& operator<<(ostream& out, const myForce& _f) {
			out << "ID :" << _f.ID << ":" << _f.name << " Const 1 : " << _f.constVal << " Const 2 : " << _f.constVal2 << " Const Vec : " << _f.constVec << " type : " << ForceType2str[_f.ftype] << endl;
			return out;
		}

		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	public:
		static unsigned int ID_gen;
		int ID;
		string name;
		double constVal;				//multiplicative constant to be applied to mass to find force
		double constVal2;
		Eigen::Vector3d constVec;				//vector constant quantity, for use with gravity
		ForceType ftype;
	};
}//namespace particleSystem
#endif