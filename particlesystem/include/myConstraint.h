#ifndef myConstraint_h
#define myConstraint_h

#include <vector>
#include <string>
#include <Eigen/Dense>
#include "myParticle.h"
#include "myForce.h"

using namespace std;

namespace particleSystem{
	class myConstraint {
	public:

		myConstraint() :ID(++ID_gen), name(""), p1(nullptr), p2(nullptr), c_Dist(), c_DistSq(), ks(), kd(), c_Type(C_NONE), useAnchor(true), anchorPoint(), p1ID(), p2ID(), p1Idx(), p2Idx(), drawCnstrPath(true) {}
		myConstraint(string _name, double _h, ConstraintType _c_Type) :ID(++ID_gen), name(_name), p1(nullptr), p2(nullptr), c_Dist(_h), c_DistSq(_h*_h),ks(), kd(), c_Type(_c_Type), useAnchor(true), anchorPoint(), p1ID(), p2ID(), p1Idx(), p2Idx(), drawCnstrPath(true) {}
		myConstraint(string _name, double _h, double _ks, double _kd, ConstraintType _c_Type, const Eigen::Ref<const Eigen::Vector3d>& _anchor) :ID(++ID_gen), name(_name), p1(nullptr), p2(nullptr), c_Dist(_h), c_DistSq(_h*_h), ks(_ks), kd(_kd), c_Type(_c_Type), useAnchor(true), anchorPoint(_anchor), p1ID(), p2ID(), p1Idx(), p2Idx(), drawCnstrPath(true) {}

		~myConstraint() {}
		//return cValue evaluated for the 2 particle positions of this constraint
		inline double calcCVal(std::shared_ptr<myParticle> p1, std::shared_ptr<myParticle> p2) {
			Eigen::Vector3d p2Pos(((p2ID != p1ID) ? p2->getPosition() : anchorPoint));			//dummy particle for anchor point	
			Eigen::Vector3d p1PosRel(p1->getPosition() - p2Pos);
			return (p1PosRel.dot(p1PosRel)) - c_DistSq; //.5 * (p1PosRel.dot(p1PosRel)) - .5 * c_DistSq;
		}
		//return cdot value evaluated for the 2 positions of this constraint
		inline double calcCDotVal(std::shared_ptr<myParticle> p1, std::shared_ptr<myParticle> p2) {
			Eigen::Vector3d p2Pos(((p2ID != p1ID) ? p2->getPosition() : anchorPoint));			//dummy particle for anchor point				
			Eigen::Vector3d p2Vel(((p2ID != p1ID) ? p2->getVelocity() : Eigen::Vector3d(0, 0, 0)));
			return (p1->getPosition() - p2Pos).dot(p1->getVelocity()  - p2Vel);									//cdot = part.pos dot part.vel
		}

		//satisfy a constraint by moving particles explicitly - need to iterate through this a few times every time step, still gonna be way faster than calculating the force constraint, and more stable
		void satisfyPosConstraint() {
			Eigen::Vector3d p1p2 = p2->getPosition() - p1->getPosition();
			double curDist = (p1p2).norm();
			Eigen::Vector3d halfCrctVec = 0.5*(p1p2*(1 - c_Dist / curDist));
			p1->modPosition(halfCrctVec);
			p2->modPosition(-1.0f * halfCrctVec);
		}

		//returns either the anchor point or the position of p1 - p1 is the "ankle" for p2->p1 constraint if no anchored
		inline Eigen::Vector3d getCurAnklePoint() { return ((p2ID != p1ID) ? p1->getPosition() : anchorPoint); }
		inline Eigen::Vector3d getInitAnklePoint() { return ((p2ID != p1ID) ? p1->initPos : anchorPoint); }


		//returns distance from initial position for "pendulum" point
		inline Eigen::Vector3d getPendDistFromInit() {
			if (p2ID != p1ID) {
				return ((p2->getPosition() - p1->getPosition()) - (p2->initPos - p1->initPos));
			} else {
				return (p2->getPosition() - p2->initPos);
			}
		}
		//returns either the anchor point or the position of p1 - p1 is the "ankle" for p2->p1 constraint if no anchored
		inline Eigen::Vector3d getCurPendVel() { if (p2ID == p1ID) { return p2->getVelocity(); }return (p2->getVelocity() - p1->getVelocity()); }
		

		//partials of c w/respect to x,y,z
		inline Eigen::Vector3d calcPartialCP1() { return p1->getPosition() - ((p2ID != p1ID) ? p2->getPosition() : anchorPoint); }			//calculate the partial derivative of C w/respect to x,y,z between 2 particles - eq; -1:id of dummy part
		inline Eigen::Vector3d calcPartialCP2() { return ((p2ID != p1ID) ? p2->getPosition() : anchorPoint) - p1->getPosition(); }														//calculate the partial derivative of C w/respect to x - eq = -partCP1Xval
		inline Eigen::Vector3d calcPartialCdotP1() { if (p2ID == p1ID) { return p1->getVelocity(); }return (p1->getVelocity() - p2->getVelocity()); }//calculate the partial derivative of C w/respect to x - eq
		inline Eigen::Vector3d calcPartialCdotP2() { if (p2ID == p1ID) { return -1*p1->getVelocity(); }return (p2->getVelocity() - p1->getVelocity()); }			//calculate the partial derivative of C w/respect to x - eq = -partCdotP1Xval

		inline Eigen::Vector3d calcPartialCP1(std::shared_ptr<myParticle> p1, std::shared_ptr<myParticle> p2) { return p1->getPosition() - ((p2ID != p1ID) ? p2->getPosition() : anchorPoint); }			//calculate the partial derivative of C w/respect to x,y,z between 2 particles - eq; -1:id of dummy part
		inline Eigen::Vector3d calcPartialCP2(std::shared_ptr<myParticle> p1, std::shared_ptr<myParticle> p2) { return -1 * calcPartialCP1(p1, p2); }														//calculate the partial derivative of C w/respect to x - eq = -partCP1Xval																																																			//partials of cdot w/respect to x,y,z
		inline Eigen::Vector3d calcPartialCdotP1(std::shared_ptr<myParticle> p1, std::shared_ptr<myParticle> p2) {if (p2ID == p1ID) { return p1->getVelocity(); }return (p1->getVelocity() - p2->getVelocity());}//calculate the partial derivative of C w/respect to x - eq
		inline Eigen::Vector3d calcPartialCdotP2(std::shared_ptr<myParticle> p1, std::shared_ptr<myParticle> p2) { return -1 * calcPartialCdotP1(p1, p2); }			//calculate the partial derivative of C w/respect to x - eq = -partCdotP1Xval

		void setP1(std::shared_ptr<myParticle> _p, int _pIdx1);
		void setP2(std::shared_ptr<myParticle> _p, int _pIdx2);



		friend ostream& operator<<(ostream& out, const myConstraint& c) {
			out << "Cnstrnt ID :" << c.ID << " Name : " << c.name << " Type : " << ConstraintType2str[c.c_Type] << " Part 1 : " << c.p1ID << " | Part 2 : " << c.p2ID;
			cout << " Anchor : " << c.anchorPoint;
			if (c.c_Type == C_Circular) { out << "\t\tConstraint radius : " << c.c_Dist; }
			out << endl;
			return out;
		}//op<<
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	public:

		static unsigned int ID_gen;
		int ID;
		string name;
		double c_Dist, c_DistSq, ks, kd;				//constraint length/radius of circle.  treat circle as constraint bar to stationary particle at center - need to choose so we can pick sq of dependent var : 
		ConstraintType c_Type;
		bool useAnchor;				//if not use anchor, uses 2 points, and draw a line between them
		Eigen::Vector3d anchorPoint;			//anchorpoint, if constraint not between 2 particles
		Eigen::Vector3d p1SPos, p1SVel, p1SAcc, p1SPosN1, p1SVelN1, p2SPos, p2SVel, p2SAcc, p2SPosN1, p2SVelN1;

		std::shared_ptr<myParticle> p1, p2;
		int p1ID, p2ID;				//particle id
		int p1Idx, p2Idx;			//index in system holder
		bool drawCnstrPath;         //if circular constraint, should i draw the path (i.e. the circle)?

	};
}//namespace particleSystem
#endif
