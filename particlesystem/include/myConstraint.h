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

		myConstraint() :ID(++ID_gen), name(""), p1(nullptr), p2(nullptr), c_Dist(), ks(), kd(), c_Type(C_NONE), useAnchor(true), anchorPoint(), p1ID(), p2ID(), p1Idx(), p2Idx(), drawCnstrPath(true) {}
		myConstraint(string _name, double _h, ConstraintType _c_Type) :ID(++ID_gen), name(_name), p1(nullptr), p2(nullptr), c_Dist(_h), ks(), kd(), c_Type(_c_Type), useAnchor(true), anchorPoint(), p1ID(), p2ID(), p1Idx(), p2Idx(), drawCnstrPath(true) {}
		myConstraint(string _name, double _h, double _ks, double _kd, ConstraintType _c_Type, const Eigen::Vector3d& _anchor) :ID(++ID_gen), name(_name), p1(nullptr), p2(nullptr), c_Dist(_h), ks(_ks), kd(_kd), c_Type(_c_Type), useAnchor(true), anchorPoint(_anchor), p1ID(), p2ID(), p1Idx(), p2Idx(), drawCnstrPath(true) {}

		~myConstraint() {}
		//return cValue evaluated for the 2 particle positions of this constraint
		inline double calcCVal(std::shared_ptr<myParticle> p1, std::shared_ptr<myParticle> p2) {
			double result;
			//Eigen::Vector3d p1Pos(p1->position[0]);
			//cout<<"p1 pos : "<<p1Pos<<endl;
			Eigen::Vector3d p2Pos(((p2ID != p1ID) ? p2->position[0] : anchorPoint));			//dummy particle for anchor point					
			result = .5 * ((p1->position[0] - p2Pos).dot(p1->position[0] - p2Pos)) - .5 * (c_Dist * c_Dist);		//circular constraint eq : c = .5 * (part.pos . part.pos) - .5 (cnstrLen * cnstrLen) - circle eq
			return result;
		}
		//return cdot value evaluated for the 2 positions of this constraint
		inline double calcCDotVal(std::shared_ptr<myParticle> p1, std::shared_ptr<myParticle> p2) {
			double result;
			//Eigen::Vector3d p1Pos(p1->position[0]);
			//Eigen::Vector3d p1Vel(p1->velocity[0]);
			Eigen::Vector3d p2Pos(((p2ID != p1ID) ? p2->position[0] : anchorPoint));			//dummy particle for anchor point				
			Eigen::Vector3d p2Vel(((p2ID != p1ID) ? p2->velocity[0] : Eigen::Vector3d(0, 0, 0)));

			result = (p1->position[0] - p2Pos).dot(p1->velocity[0] - p2Vel);									//cdot = part.pos dot part.vel
			return result;
		}
		////return cValue evaluated for the 2 particle positions of this constraint-spd implementation
		//inline double calcCValSPD(std::shared_ptr<myParticle> p1, std::shared_ptr<myParticle> p2, double delT) {
		//	double result;
		//	p1SPos = p1->position[0];
		//	p1SVel = p1->velocity[0];
		//	//cout<<"p1 pos : "<<p1Pos<<endl;
		//	p2SPos = ((p2ID != p1ID) ? p2->position[0] : anchorPoint);
		//	p2SVel = ((p2ID != p1ID) ? p2->velocity[0] : Eigen::Vector3d(0, 0, 0));

		//	p1SPosN1 = p1SPos + delT*p1SVel;
		//	p2SPosN1 = p2SPos + delT*p2SVel;

		//	result = .5 * ((p1SPosN1 - p2SPosN1).dot(p1SPosN1 - p2SPosN1)) - .5 * (c_Dist * c_Dist);		//circular constraint eq : c = .5 * |(part.pos . part.pos)| - .5 (cnstrLen * cnstrLen) <-- circle eq
		//	//result = (p1PosN1 - p2PosN2;		//circular constraint eq : c = .5 * |(part.pos . part.pos)| - .5 (cnstrLen * cnstrLen) <-- circle eq
		//	return result;
		//}
		////return cdot value evaluated for the 2 positions of this constraint - spd implementation
		//inline double calcCDotValSPD(std::shared_ptr<myParticle> p1, std::shared_ptr<myParticle> p2, double delT) {
		//	double result;
		//	p1SPos = p1->position[0];
		//	p1SVel = p1->velocity[0];
		//	p1SAcc = p1->forceAcc[0];
		//	p2SPos = (p2ID != p1ID ? p2->position[0] : anchorPoint);			//dummy particle for anchor point				
		//	p2SVel = (p2ID != p1ID ? p2->velocity[0] : Eigen::Vector3d(0, 0, 0));
		//	p2SAcc = (p2ID != p1ID ? p2->forceAcc[0] : Eigen::Vector3d(0, 0, 0));

		//	Eigen::Vector3d p1SPosN1(p1SPos + delT*p1SVel), p2SPosN1(p2SPos + delT*p2SVel);
		//	////d vec from SPD literature
		//	//Eigen::Vector3d dVec(p1SPosN1 - p2SPosN1);
		//	//dVec.normalize();
		//	////D matrix from literature
		//	//Mat3d ddT = oprod(dVec,dVec);

		//	Eigen::Vector3d p1VelN1(p1SVel + delT*p1SAcc), p2VelN1(p2SVel + delT*p2SAcc);

		//	result = ((p1SPosN1 - p2SPosN1).dot(p1VelN1 - p2VelN1));									//cdot = part.pos dot part.vel
		//	return ((p1SPosN1 - p2SPosN1).dot(p1VelN1 - p2VelN1));
		//}
		//satisfy a constraint by moving particles explicitly - need to iterate through this a few times every time step, still gonna be way faster than calculating the force constraint, and more stable
		void satisfyPosConstraint() {
			Eigen::Vector3d p1p2 = p2->position[0] - p1->position[0];
			double curDist = (p1p2).norm();
			Eigen::Vector3d halfCrctVec = 0.5*(p1p2*(1 - c_Dist / curDist));
			p1->position[0] += halfCrctVec;
			p2->position[0] -= halfCrctVec;
		}

		//partials of c w/respect to x,y,z
		inline Eigen::Vector3d calcPartialCP1(std::shared_ptr<myParticle> p1, std::shared_ptr<myParticle> p2) { return p1->position[0] - ((p2ID != p1ID) ? p2->position[0] : anchorPoint); }			//calculate the partial derivative of C w/respect to x,y,z between 2 particles - eq; -1:id of dummy part
		inline Eigen::Vector3d calcPartialCP2(std::shared_ptr<myParticle> p1, std::shared_ptr<myParticle> p2) { return -1 * calcPartialCP1(p1, p2); }														//calculate the partial derivative of C w/respect to x - eq = -partCP1Xval
		//partials of cdot w/respect to x,y,z
		inline Eigen::Vector3d calcPartialCdotP1(std::shared_ptr<myParticle> p1, std::shared_ptr<myParticle> p2) {
			//Eigen::Vector3d p1Pos(p1->position[0]);
			//Eigen::Vector3d p1Vel(p1->velocity[0]);
			//Eigen::Vector3d p1Force(p1->forceAcc[0]);
			//Eigen::Vector3d p2Pos( ((p2->getID() != p1->getID() ) ? p2->position[0] : anchorPoint));			//dummy particle for anchor point				
			Eigen::Vector3d p2Vel(((p2ID != p1ID) ? p2->velocity[0] : Eigen::Vector3d(0, 0, 0)));

			return (p1->velocity[0] - p2Vel);					//cdot = (part.pos.x * part.vel.x) + (part.pos.y * part.vel.y) + (part.pos.z * part.vel.z)
		}//calculate the partial derivative of C w/respect to x - eq

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
		double c_Dist, ks, kd;				//constraint length/radius of circle.  treat circle as constraint bar to stationary particle at center - need to choose so we can pick sq of dependent var : 
		ConstraintType c_Type;
		bool useAnchor;				//if not use anchor, uses 2 points, and draw a line between them
		Eigen::Vector3d anchorPoint;			//anchorpoint, if constraint not between 2 particles
		Eigen::Vector3d p1SPos, p1SVel, p1SAcc, p1SPosN1, p1SVelN1, p2SPos, p2SVel, p2SAcc, p2SPosN1, p2SVelN1;
		//myParticle *p1, *p2;	
		std::shared_ptr<myParticle> p1, p2;
		int p1ID, p2ID;				//particle id
		int p1Idx, p2Idx;			//index in system holder
		bool drawCnstrPath;         //if circular constraint, should i draw the path (i.e. the circle)?

	};
}//namespace particleSystem
#endif
