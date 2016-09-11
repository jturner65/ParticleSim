#ifndef myCollider_h
#define myCollider_h

#include <vector>
#include <string>
#include <Eigen/Dense>
#include "myParticle.h"
#include <memory>

using namespace std;

namespace particleSystem{

	class myCollider {
	public:
		myCollider() :ID(++ID_gen), name(), colType(), drawLoc(), center(), radius(), minMaxRadius(2), planeNormal(), verts(), peq(4), intRefl(), Krest(1), muFrict(0) {}

		myCollider(string _n, const Eigen::Ref<const Eigen::Vector3d>& _drawLoc, const Eigen::Ref<const Eigen::Vector3d>& _ctr, const Eigen::Ref<const Eigen::Vector3d>& _rad, bool _inRefl) :							//sphere collider
			ID(++ID_gen), name(_n), colType(SPHERE), drawLoc(_drawLoc), center(_ctr), radius(_rad), minMaxRadius(2), planeNormal(), verts(), peq(4), intRefl(_inRefl), Krest(1) {
			initCollider();
		}

		myCollider(string _n, const Eigen::Ref<const Eigen::Vector3d>& _drawLoc, vector<Eigen::Vector3d> _vs, bool _inRefl) :		//flat collider
			ID(++ID_gen), name(_n), colType(FLAT), drawLoc(_drawLoc), center(), radius(), minMaxRadius(2), planeNormal(), verts(_vs), peq(4), intRefl(_inRefl), Krest(1) {
			initCollider();
		}

		~myCollider() {};

		int checkCollision(double deltaT, std::shared_ptr<myParticle> part);		//0 vec if no collision, 1 if collision via breach, 2 if collision next timestep, 3 if need to counter force due to contact

		void initCollider();

		Eigen::Vector3d handlePlanarCollision(std::shared_ptr<myParticle> p, int res);
		void handleSphereCollision(std::shared_ptr<myParticle> p, int res);
		void buildPlaneNorm();		//build normal for this flat object given verts exist
		void findPlaneEQ();
		void findMinMaxRadius();
		vector<Eigen::Vector3d> getPartVelNorm(const Eigen::Ref<const Eigen::Vector3d>& partVel, const Eigen::Ref<const Eigen::Vector3d>& norm);

		void handleVerletCol(std::shared_ptr<myParticle> p);

		inline void setNormal(const Eigen::Ref<const Eigen::Vector3d>& _n) { planeNormal = _n; planeNormal.normalize(); }

		inline double getDistFromCenter(const Eigen::Ref<const Eigen::Vector3d>& _p, double mult) { return ((_p - center) * (mult == 0 ? 1 : ((_p[2] < center[2]) ? -mult : mult))).norm(); }

		//inline Eigen::Vector3d getNormal(){	return planeNormal;}//getNormal for flat plane - need to have verts and precalculated normal

		inline Eigen::Vector3d getSphereNormal(const Eigen::Ref<const Eigen::Vector3d>& _loc) {//get normal at a particular location - no matter where inside or outside of sphere, normal built from this point and center will point in appropriate dir
			//if sphere, normal will be either pointing out or in, colinear with line from center to _loc
			Eigen::Vector3d normDir = (center - _loc);//either point into center if internal reflections or point out of center if not
			double mult = ((intRefl) ? 1 : -1);
			normDir *= mult;
			normDir.normalize();
			return normDir;
		}//getNormal
		inline bool isIntRefl() { return intRefl; }

		inline double min3(double x, double y, double z) { return (x < y ? (x < z ? x : z) : (y < z ? y : z)); }
		inline double max3(double x, double y, double z) { return (x > y ? (x > z ? x : z) : (y > z ? y : z)); }

		friend ostream& operator<<(ostream& out, const myCollider& col) {
			out << "Col ID: " << col.ID << " Name : " << col.name << " Type : " << CollisionType2str[col.colType];
			if (col.colType == FLAT) {
				out << " Normal " << col.planeNormal << " num verts : " << col.verts.size();
				for (unsigned int i = 0; i < col.verts.size(); ++i) { out << "\n\tvert[" << i << "]:" << col.verts[i] << " | "; }
				for (unsigned int i = 0; i < 4; ++i) { out << "\n\tpeq[" << i << "]:" << col.peq[i] << " | "; }
			}
			else if (col.colType == SPHERE) {
				out << " Center : " << col.center << " Radius : " << col.radius << " internal relfections?" << col.intRefl;
			}
			out << "Krest : " << col.Krest << endl;
			return out;
		}//op <<
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	public:
		static unsigned int ID_gen;
		int ID;
		string name;
		CollisionType colType;
		Eigen::Vector3d drawLoc;			//drawn location of this collider - only different from center if needed for display purposes
		//shere
		Eigen::Vector3d center;			//actual location of this collider, w/respect to particles - used for sphere
		Eigen::Vector3d radius;			//radius around center for ellipsoid (in each direction for x,y,z)
		vector<double> minMaxRadius;		//minimum dist from center to still be inside (in case of ellipse, to minimize calculations)
		//plane
		Eigen::Vector3d planeNormal;		//normal of this collider, if flat plane
		vector<Eigen::Vector3d> verts;	//vertices of this object, if flat plane
		vector<double> peq;		//plane equation values
		bool  intRefl;			//internal reflections? for sphere (collide on inside)
		double Krest;			//coefficent of restitution - how much bounce do we have :1 total bounce, 0 no bounce.  multiply this against part's Velperp
		double muFrict;         //friction coefficient
	};//myCollider
}//namespace particleSystem

#endif
