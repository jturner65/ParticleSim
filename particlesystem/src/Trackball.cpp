//////////////////////////////////////////////////////////////////////
// 
//	File: Trackball.cpp
// 
//	Graphics and Imaging Lab
//	Dept. of Computer Scinece and Engineering
//	University of Washington
//
// 
//	Created:       Thu Aug 04 13:48 2001 by C. Karen Liu
//	
//	Trackball definitions
//////////////////////////////////////////////////////////////////////

#ifndef	__TRACKBALL_H__
#include "Trackball.h"
#endif	//__TRACKBALL_H__


Eigen::Vector3d Trackball::MouseOnSphere(double mouseX, double mouseY)
{   
	double mag;
	Eigen::Vector3d pointOnSphere;

    pointOnSphere[0] = (mouseX - mCenter[0]) / mRadius;
    pointOnSphere[1] = (mouseY - mCenter[1]) / mRadius;
    
	mag = pointOnSphere[0] * pointOnSphere[0] + pointOnSphere[1]*pointOnSphere[1];
    if (mag > 1.0) 
	{	register double scale = 1.0/sqrt(mag);
		pointOnSphere[0] *= scale; 
		pointOnSphere[1] *= scale;
		pointOnSphere[2] = 0.0;
    } 
	else 
		pointOnSphere[2] = sqrt(1 - mag);
    
	return pointOnSphere;
}

Eigen::Vector4d Trackball::Qt_Conj(Eigen::Vector4d& q)
{
	Eigen::Vector4d qq;
    qq[0] = -q[0]; qq[1] = -q[1]; qq[2] = -q[2]; qq[3] = q[3];
    return (qq);
}


Eigen::Vector4d Trackball::Qt_FromBallPoints(Eigen::Vector3d& from, Eigen::Vector3d& to)
{
	Eigen::Vector4d tempVec4;
    tempVec4[0] = from[1]*to[2] - from[2]*to[1];
    tempVec4[1] = from[2]*to[0] - from[0]*to[2];
    tempVec4[2] = from[0]*to[1] - from[1]*to[0];
    tempVec4[3] = from[0]*to[0] + from[1]*to[1] + from[2]*to[2];
	return tempVec4;
}


Eigen::Vector4d Trackball::Qt_Mul(Eigen::Vector4d& qL, Eigen::Vector4d& qR)
{
	Eigen::Vector4d qq;
    qq[3] = qL[3]*qR[3] - qL[0]*qR[0] - qL[1]*qR[1] - qL[2]*qR[2];
    qq[0] = qL[3]*qR[0] + qL[0]*qR[3] + qL[1]*qR[2] - qL[2]*qR[1];
    qq[1] = qL[3]*qR[1] + qL[1]*qR[3] + qL[2]*qR[0] - qL[0]*qR[2];
    qq[2] = qL[3]*qR[2] + qL[2]*qR[3] + qL[0]*qR[1] - qL[1]*qR[0];
    return qq;
}

void Trackball::Qt_ToMatrix(Eigen::Vector4d& q, Eigen::Matrix4d &out)
{
    double Nq = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    double s = (Nq > 0.0) ? (2.0 / Nq) : 0.0;
    double xs = q[0]*s,	      ys = q[1]*s,	  zs = q[2]*s;
    double wx = q[3]*xs,	      wy = q[3]*ys,	  wz = q[3]*zs;
    double xx = q[0]*xs,	      xy = q[0]*ys,	  xz = q[0]*zs;
    double yy = q[1]*ys,	      yz = q[1]*zs,	  zz = q[2]*zs;
    out(0,0) = 1.0 - (yy + zz); 
	out(1,0) = xy + wz; 
	out(2,0) = xz - wy;
    out(0,1) = xy - wz; 
	out(1,1) = 1.0 - (xx + zz); 
	out(2,1) = yz + wx;
    out(0,2) = xz + wy; 
	out(1,2) = yz - wx; 
	out(2,2) = 1.0 - (xx + yy);
    out(0,3) = out(1,3) = out(2,3) = out(3,0) = out(3,1) = out(3,2) = 0.0;
    out(3,3) = 1.0;
}

void Trackball::StartBall(double x, double y)
{
	mStartPix = Eigen::Vector2d(x, y);
	mStartPos = MouseOnSphere(x,y);
}


void Trackball::UpdateBall(double x, double y)
{
	if(mNonCenter) {
		x -= mStartPix[0];
		y -= mStartPix[1];
	}
	Eigen::Vector3d toPos = MouseOnSphere(x,y);
	Eigen::Vector4d newQuat = Qt_FromBallPoints(mStartPos, toPos);
	mStartPos = toPos;
	mCurrQuat = Qt_Mul(newQuat, mCurrQuat);
}


void Trackball::GetCurrentRotation(double *returnRot)
{	
	Eigen::Matrix4d tempRot;
	Qt_ToMatrix(mCurrQuat,tempRot);
	for(int x=0;x<4;x++)
	{	for(int y=0;y<4;y++)
		{returnRot[4*x+y]=tempRot(y,x);}
	}
}

void Trackball::GetCurrentRotation(Eigen::Matrix4d & returnRot)
{	
	Qt_ToMatrix(mCurrQuat,returnRot);

}

void Trackball::GetCurrentRotation(Eigen::Vector4d & returnQuat) {
	returnQuat = mCurrQuat;
}

void Trackball::Draw(int winWidth, int winHeight)
{
	// draw the axes
	glPushMatrix();
	if(mNonCenter) glTranslated(mCenter[0], mCenter[1], mCenter[2]);
	double *rotationMat = new double[16];
	GetCurrentRotation(rotationMat);
	glMultMatrixd( rotationMat );
	glDisable( GL_LIGHTING );
		glLineWidth(2.0);
		glColor3f( 1.0f, 0.0f, 0.0f );
		double len = 0.2;
		glBegin( GL_LINES );
		glVertex3d( -0.0f, 0.0f, -0.0f );
		glVertex3d( len, 0.0f, -0.0f );
		glEnd();

		glColor3f( 0.0f, 1.0f, 0.0f );
		glBegin( GL_LINES );
		glVertex3d( 0.0f, -0.0f, 0.0f );
		glVertex3d( 0.0f, len, 0.0f );
		glEnd();

		glColor3f( 0.0f, 0.0f, 1.0f );
		glBegin( GL_LINES );
		glVertex3d( 0.0f, 0.0f, -0.0f );
		glVertex3d( 0.0f, 0.0f, len );
		glEnd();
	glEnable( GL_LIGHTING );
	glPopMatrix();

	// draw the globe
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glViewport( 0, 0, winWidth, winHeight );
	gluOrtho2D(0.0, (GLdouble)winWidth, 0.0, (GLdouble)winHeight);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glEnable(GL_BLEND);
	glShadeModel(GL_SMOOTH);
	
    glDisable( GL_LIGHTING );
	glLineWidth(2.0);
	if(mNonCenter) glTranslated(mStartPix[0], mStartPix[1], 0);
	glColor3f( 1.0f, 1.0f, 0.0f );
		double r = mRadius;
		glBegin( GL_LINE_LOOP);		
		for(int i = 0; i < 360; i+=4){
			double theta = i / 180.0 * 3.14156;
			double x = r * cos(theta);
			double y = r * sin(theta);
			glVertex2d( (GLdouble)((winWidth >> 1) + x), (GLdouble)((winHeight >> 1) + y));
		}
		glEnd();
   
	glPopMatrix();
    glEnable( GL_LIGHTING );
}
