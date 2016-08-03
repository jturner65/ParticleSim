#ifndef mySpring_h
#define mySpring_h

#include <vector>
#include <string>
#include <iostream>
#include <Eigen/Dense>

#include "myParticle.h"
#include "myForce.h"

using namespace std;

namespace particleSystem{

	class mySpring {
	public:
		mySpring(void) :ID(++ID_gen), name(""), Ks(0), Kd(0), restLen(0), a(nullptr), b(nullptr), Jp(), Jv(), Id3x3() { Jp.setZero(); Id3x3.setIdentity(); }
		mySpring(string _n, std::shared_ptr<myParticle> _a, std::shared_ptr<myParticle> _b, double _k1, double _k2, double _rstl) :ID(++ID_gen), name(_n), Ks(_k1), Kd(_k2), restLen(_rstl), baseRestLen(_rstl), a(_a), b(_b), Jp(), Jv(), Id3x3() { Jp.setZero(); Id3x3.setIdentity(); initSpring(); }
		mySpring(const mySpring& other) :ID(other.ID), name(other.name), Ks(other.Ks), Kd(other.Kd), restLen(other.restLen), a(other.a), b(other.b), Jp(other.Jp), Jv(other.Jv), Id3x3() { Id3x3.setIdentity(); }
		~mySpring(void) {}

		inline void initSpring() {
			Jv.setIdentity();       //ident
			Jv *= Kd;
			buildSprJpJv();
		}
		inline void scaleSpring(double mult) {
			Ks *= mult;
			Kd *= mult;
			initSpring();
		}

		inline void calcAndApplySpringForce() {
			vector<Eigen::Vector3d> result = calcSpringForce();
			a->applyForce(result[0]);
			b->applyForce(result[1]);
		}

		void buildSprJpJv();
		void reCalcRestLen() { restLen = (a->position[0] - b->position[0]).norm(); }
		vector<Eigen::Vector3d> calcSpringForce();//sprRestVec is vector in direction from p2 to p1 at rest state (rest length)

		friend ostream& operator<<(ostream& out, const mySpring& _f) {
			out << "ID :" << _f.ID << ":" << _f.name << " Ks : " << _f.Ks << " Kd : " << _f.Kd << endl;
			return out;
		}
		void setRestLenByMult(double mult) {
			restLen = mult*baseRestLen + baseRestLen;
		}

	public:
		static unsigned int ID_gen;
		int ID;
		string name;
		std::shared_ptr<myParticle> a, b;       //particles on either end of this spring

		Eigen::Matrix3d Jp, Jv, Id3x3;                           //spring jacobians for implicit calculation

		double Ks;				                //multiplicative constant to be applied to mass to find force
		double Kd;
		double restLen, baseRestLen;                         //rest length of this spring

	};
}//namespace particleSystem
#endif