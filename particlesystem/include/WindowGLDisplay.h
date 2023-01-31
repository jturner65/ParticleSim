#ifndef _WINDOWGLDISPLAY_
#define _WINDOWGLDISPLAY_

#include "Trackball.h"

#include <fltk/run.h>
#include <fltk/GlWindow.h>
#include <fltk/draw.h>
#include <fltk/File_Chooser.h>
#include <fltk/gl.h>
#include <fltk/glut.h>
#ifdef WIN32
#include <windows.h>
#include <GL/glu.h>
#elif CYGWIN
#include <GL/glu.h>
#elif MACOSX
#include <opengl/glu.h>
#endif
#include <stdlib.h>

#include <Eigen/Dense>

#include <vector>
#include <memory>

namespace particleSystem{
    class myParticle;
    class myForce;
    class mySystem;
    class myConstraint;
    class mySpring;
    class myCollider;
}

using namespace particleSystem;

class WindowGLDisplay : public fltk::GlWindow
{
public:
	WindowGLDisplay(int x, int y, int w, int h);

	bool mAntiAlias;
	//virtual functions that replace GlWindow
	void draw();
	int FindPartIDXByID(std::vector<std::shared_ptr<myParticle>>& partAra, int id);
	int FindCnstrntIDXByID(std::vector<std::shared_ptr<myConstraint>>& cnstrntAra, int id);
	void drawFluidVel(std::shared_ptr<mySystem> system);
	void drawParts(std::vector<std::shared_ptr<particleSystem::myParticle>>& partAra, double calc_partSize, bool draw1stPartBlue);
	void drawPartsCOM(double calc_partSize);
	//void drawCnstrntLine(std::vector<std::shared_ptr<particleSystem::myParticle>>& partAra, std::vector<std::shared_ptr<particleSystem::myConstraint>>& cnstrntAra, std::vector<std::vector<int>>& linePointIDS, bool bindCnstr, double cnstrR);	
	void drawCnstrnt(std::vector<std::shared_ptr<particleSystem::myParticle>>& partAra, std::vector<std::shared_ptr<particleSystem::myConstraint>>& cnstrntAra);
	void drawSpring(std::vector<std::shared_ptr<particleSystem::mySpring>>& springAra);
	void render();
	int handle(int event);
	void refresh();
	void resizeWindow(int width, int height);
	void resize( int x, int y, int w, int h );

	bool mShowConstraints;
	bool mShowModel;
	bool mShowMarker;
	bool mPlayback;

	//precalced sin and cos arrays
	float cosAra[629], sinAra[629];


	int mWinWidth, mWinHeight;
	Trackball mTrackball;
	bool mMouseDown;
	// void* mContext;
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
	void resetCamera();

	int mSelectedMark;
	Eigen::Vector3d mMouse;
	double mXstart, mYstart, mXangle, mYangle, mZangle, mZoom, mTranslateX, mTranslateY;
	bool mControlState, mShiftState, mAltState, mTabState, mLockState;
	Eigen::Vector3d mFixedTrans;
	bool mRotatingXY, mZooming, mTranslating, mPicking;
	std::string mCommand;
	int mWhichClick;
	Eigen::Vector3d mForceDir;

	Eigen::Vector3d mInitEye;
	Eigen::Vector3d mInitLookAt;
	Eigen::Vector3d mInitUp;

	double mScaleGL;	// pixel to coord frame scale for 2D

	// picking stuff
	Trackball *mBodyTrackBall;

	// gets new 3D position of point in same plane as point "pos"
	Eigen::Vector3d getNew3DPos(const Eigen::Ref<const Eigen::Vector3d>& pos);
};

#endif
