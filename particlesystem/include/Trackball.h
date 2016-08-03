//////////////////////////////////////////////////////////////////////
// 
//	File: Trackball.h
// 
//	Graphics and Imaging Lab
//	Dept. of Computer Scinece and Engineering
//	University of Washington
//
// 
//	Created:       Thu Aug 04 13:47 2001 by C. Karen Liu
//	
//	Trackball controls the user interaction 
//////////////////////////////////////////////////////////////////////

#ifndef __TRACKBALL_H__
#define __TRACKBALL_H__


#include <Eigen/Dense>

#ifdef WIN32
#include <windows.h>
#include <GL/glu.h>
#elif CYGWIN
#include <GL/glu.h>
#elif MACOSX
#include <opengl/glu.h>
#endif

//Trackball Routine

class Trackball
{
	public:
		Trackball(Eigen::Vector3d& center, double radius){mCenter = center;mRadius=radius; mCurrQuat = Eigen::Vector4d(0, 0, 0, 1);mNonCenter=false;}
		void Ball_Place(Eigen::Vector3d& center, double radius){mCenter = center;mRadius=radius;}
		void Ball_ComputeRotationMatrix(Eigen::Matrix4d &);
		void StartBall(double x, double y);
		void UpdateBall(double x, double y);
		void GetCurrentRotation(double *returnRot);
		void GetCurrentRotation(Eigen::Matrix4d & returnRot);
		void GetCurrentRotation(Eigen::Vector4d& returnQuat);
		void Draw(int, int);
		void SetCenter(Eigen::Vector3d& cent){ mCenter = cent;}
		void SetCurrQuat(Eigen::Vector4d& q){ mCurrQuat = q; }
		double mRadius;
		inline Eigen::Vector3d& getCenter(){return mCenter;}
		inline void setNoCenter(){mNonCenter=true;}

		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	private:	
		void Qt_ToMatrix(Eigen::Vector4d& q, Eigen::Matrix4d &out);
		Eigen::Vector3d MouseOnSphere(double mouseX, double mouseY);
		Eigen::Vector4d Qt_FromBallPoints(Eigen::Vector3d& from, Eigen::Vector3d& to);
		Eigen::Vector4d Qt_Mul(Eigen::Vector4d& qL, Eigen::Vector4d& qR);
		Eigen::Vector4d Qt_Conj(Eigen::Vector4d& q);
		
		Eigen::Vector2d mStartPix;
		bool mNonCenter;

		Eigen::Vector3d mCenter;
		Eigen::Vector3d mStartPos;
		Eigen::Vector4d mCurrQuat;
};

#endif
