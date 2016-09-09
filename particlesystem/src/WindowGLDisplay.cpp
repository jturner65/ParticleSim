#include "ControlPanel.h"
#include "WindowGLDisplay.h"

#include <Eigen/Dense>
#include "glFuncs.h"
#include "MyParticleWorld.h"
#include "ParticleSystemUI.h"
#include "myParticle.h"
#include "myForce.h"
#include "mySolver.h"
#include "mySystem.h"
#include "myConstraint.h"
using namespace particleSystem;

#define BUFSIZE 1024
GLuint selectBuf[BUFSIZE];
GLint hits;
extern int gCurrentScene;
extern ParticleSystemUI* mUI;


void LightInit()
{

	static float ambient[]             = {0.3f, 0.3f, 0.3f, 1.0f};
	static float diffuse[]             = {0.6f, 0.6f, 0.6f, 1.0f};
	static float front_mat_shininess[] = {60.0f};
	static float front_mat_specular[]  = {0.2f, 0.2f,  0.2f,  1.0f};
	static float front_mat_diffuse[]   = {0.5f, 0.28f, 0.38f, 1.0f};
	static float lmodel_ambient[]      = {0.2f, 0.2f,  0.2f,  1.0f};
	static float lmodel_twoside[]      = {GL_FALSE};
	if ((SNOW_GLOBE == gCurrentScene) || (MPM_FLUID == gCurrentScene) ){
		ambient[0] = 0.8f; ambient[1] = 0.8f; ambient[2] = 0.8f; 
		diffuse[0] = 1.0f; diffuse[1] = 1.0f; diffuse[2] = 1.0f;

	}
	GLfloat position[] = {1.0,0.0,0.0,0.0};
	GLfloat position1[] = {-1.0,0.0,0.0,0.0};
	GLfloat position2[] = {0.0,0.0,1.0,0.0};
	GLfloat position3[] = {0.0,0.0,-1.0,0.0};

	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_AMBIENT,  ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, position);

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT,  lmodel_ambient);
	glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, lmodel_twoside);

	glEnable( GL_LIGHT1);
	glLightfv(GL_LIGHT1,GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT1,GL_POSITION, position1);
	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);

	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, front_mat_shininess);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  front_mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   front_mat_diffuse);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_CULL_FACE);
	glEnable(GL_NORMALIZE);
}

void inits(void){
	//glClearColor(242.0f / 255, 247.0f / 255, 1.0f, 1.0f);
	if (SNOW_GLOBE == gCurrentScene) {
		glClearColor(0, 0, 0, 1.0f);
	}
	else if (MPM_FLUID == gCurrentScene) {
		glClearColor(0, 0, .2f, 1.0f);
	}
	else {
		glClearColor(242.0f / 255, 247.0f / 255, 1.0f, 1.0f);
	}
	glClearAccum(0.0f, 0.0, 0.0, 0.0);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glShadeModel(GL_SMOOTH);
	LightInit();
}


WindowGLDisplay::WindowGLDisplay(int x, int y, int w, int h): GlWindow(x,y,w,h,NULL),mTrackball(Eigen::Vector3d(0,0,0), h/2.5){
	mWinWidth = w;
	mWinHeight = h;
	mSelectedMark = -1;

	mZoom = mZangle = mTranslateX = mTranslateY = 0;

	mXangle = mXstart = 0;
	mYangle = mYstart = 180;
	mControlState = mShiftState = mAltState = mLockState = mTabState = false;
	mPicking = mRotatingXY = mZooming = mTranslating = false;
	mFixedTrans = Eigen::Vector3d(0,0,0);
	mMouseDown = false;

	mShowConstraints = true;
	mShowModel =  true;
	mShowMarker = false;
	mPlayback = false;

	mForceDir = Eigen::Vector3d(0, 0, 0);

	mInitEye = Eigen::Vector3d(0,0,3.0);
	mInitLookAt = Eigen::Vector3d(0, 0, -100.0) ;
	mInitUp = Eigen::Vector3d(0, 1, 0) ;

	mScaleGL = 0.012;

	mBodyTrackBall=NULL;
	mAntiAlias = false;
}

void WindowGLDisplay::draw(){
	inits();
	resetCamera();
	render();
}


void WindowGLDisplay::render()
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	double *rotationMat = new double[16];
	mTrackball.GetCurrentRotation(rotationMat);
	glMultMatrixd( rotationMat );

	glDisable( GL_LIGHTING );
	if(mRotatingXY || mZooming || mTranslating){
		glColor3f( 1.0f, 0.0f, 0.0f );
		DrawArrow(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 0), 0.4, 0.01);

		glColor3f( 0.0f, 1.0f, 0.0f );
		DrawArrow(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 0), 0.4, 0.01);

		glColor3f( 0.0f, 0.0f, 1.0f );
		DrawArrow(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1), 0.4, 0.01);
	}
	glEnable( GL_LIGHTING );

	glTranslated(mFixedTrans[0] / 30.0, mFixedTrans[1] / 30, mFixedTrans[2] / 10);			//moving camera around. YUCK bad template code bad.

	// curr frame to be displayed
	unsigned int currFrame = mUI->mControl->getCurrentFrame();

	// ------------------------------ your drawing code starts here  --------------------------------
	
	vector<std::shared_ptr<particleSystem::myParticle>> partAra;
	vector<std::shared_ptr<particleSystem::mySpring>> springAra;
	vector<std::shared_ptr<particleSystem::myConstraint>> cnstrntAra = MyParticleWorld::systemConstraintHolder[gCurrentScene][currFrame];

	int pSize = MyParticleWorld::systems[gCurrentScene]->p.size();

    //decide whether drawing for during simulation or during playback and set this to appropriate idx
	int stIdx = ((mUI->mSim_but->value()) ? MyParticleWorld::systemParticleHolder[gCurrentScene].size() - 1 : currFrame);
	partAra = MyParticleWorld::systemParticleHolder[gCurrentScene][stIdx];
	if (currFrame<MyParticleWorld::systemConstraintHolder[gCurrentScene].size()) { cnstrntAra = MyParticleWorld::systemConstraintHolder[gCurrentScene][stIdx]; }
	if (currFrame<MyParticleWorld::systemSpringHolder[gCurrentScene].size()) { springAra = MyParticleWorld::systemSpringHolder[gCurrentScene][stIdx]; }
	double partRadDraw = (partRad), gndCubeDim = 10.0;

    //draw scene-specific stuff
	if ((SNOW_GLOBE == gCurrentScene) || (MPM_FLUID == gCurrentScene)) {
		partRadDraw *= .5;
		//TODO move collider drawing to colliders.  derp
		glPushMatrix();//draw sphere
			glColor4d(0.3f, 0.7f, 0.7f, 0.1f);
			glTranslated(MyParticleWorld::systems[gCurrentScene]->colliders[1]->center(0), MyParticleWorld::systems[gCurrentScene]->colliders[1]->center(1), MyParticleWorld::systems[gCurrentScene]->colliders[1]->center(2));
			glutWireSphere(snowGlobRad, 32, 32);
		glPopMatrix();
	}
	if((BALL_DROP ==gCurrentScene) || (INV_PEND ==gCurrentScene) || (SNOW_GLOBE == gCurrentScene) || (MPM_FLUID == gCurrentScene) || (SEAWEED == gCurrentScene) || (MSPR_MTN_PROJ == gCurrentScene)|| (MSPR_MTN_PRO2 == gCurrentScene)){
		glPushMatrix();//draw ground
			glColor4d(0.3, 0.7, 0.7, 0.1);
			glTranslated(MyParticleWorld::systems[gCurrentScene]->colliders[0]->drawLoc(0), MyParticleWorld::systems[gCurrentScene]->colliders[0]->drawLoc(1) - (.5*gndCubeDim), MyParticleWorld::systems[gCurrentScene]->colliders[0]->drawLoc(2));
			glScalef(5,1,5);
			glutSolidCube(gndCubeDim); //ground
		glPopMatrix();
    }
    if((MSPR_MTN_PROJ == gCurrentScene) || (MSPR_MTN_PRO2 == gCurrentScene)){
        drawSpring(springAra);
    }
	drawCnstrnt(partAra, cnstrntAra);
	drawParts(partAra, partRadDraw * 2, partRadDraw * 2, (gCurrentScene == BALL_DROP));
	if (MyParticleWorld::systems[gCurrentScene]->flags[mySystem::showVel] && ((SNOW_GLOBE == gCurrentScene) || (MPM_FLUID == gCurrentScene) || (SEAWEED == gCurrentScene))) {
        drawFluidVel( MyParticleWorld::systems[gCurrentScene]);
    }

	// -------------------------------your drawing code ends here  -----------------------------------

	// draw the editing ball
	if(mBodyTrackBall) {
		mBodyTrackBall->Draw(mWinWidth, mWinHeight);
	}

	// Draw the globe
	if(mRotatingXY){
		glPushMatrix();
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glViewport( 0, 0, mWinWidth, mWinHeight );
		gluOrtho2D(0.0, (GLdouble)mWinWidth, 0.0, (GLdouble)mWinHeight);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glDisable( GL_LIGHTING );
		glLineWidth(2.0);
		glColor3f( 1.0f, 1.0f, 0.0f );
		double r = mTrackball.mRadius;
		glBegin( GL_LINE_LOOP);		
		for(int i = 0; i < 360; i+=4){
			double theta = i / 180.0 * M_PI;
			double x = r * cos(theta);
			double y = r * sin(theta);
			glVertex2d( (GLdouble)((mWinWidth >> 1) + x), (GLdouble)((mWinHeight >> 1) + y));
		}
		glEnd();
		glPopMatrix();
		glEnable( GL_LIGHTING );
	}
	delete [] rotationMat;

	glFlush();
}

int WindowGLDisplay::FindPartIDXByID(vector<std::shared_ptr<myParticle>>& partAra, int id){for(unsigned int i = 0; i < partAra.size(); ++i){if (partAra[i]->ID == id) return i;} return -1;}
int WindowGLDisplay::FindCnstrntIDXByID(vector<std::shared_ptr<myConstraint>>& cnstrntAra, int id){for(unsigned int i = 0; i < cnstrntAra.size(); ++i){if (cnstrntAra[i]->ID == id) return i;} return -1;}

//draw fluid vectors if a fluid box exists
void WindowGLDisplay::drawFluidVel(shared_ptr<mySystem> system){
        //draw force lines for fluid boxes        
    glPushMatrix();
    glColor4d(1.0f, 0.0f, 0.0f, .1f);
	Eigen::Vector3d v(0,0,0);
    float vecLen = 1000;
	glTranslated(system->fluidBox->startLoc(0),	 system->fluidBox->startLoc(1), system->fluidBox->startLoc(2));
	Eigen::Vector3d vertLoc(0, 0, 0);
	for (unsigned int k = 0; k < system->fluidBox->sx1i; ++k) {
		for (unsigned int j = 0; j < system->fluidBox->sy1i; ++j) {
			for (unsigned int i = 0; i < system->fluidBox->sz1i; ++i) {
				vertLoc << i + .5, j + .5, k + .5;
				vertLoc.cwiseProduct(system->fluidBox->cellSz);
				glPushMatrix();
				glTranslated(vertLoc(0), vertLoc(1), vertLoc(2));

				glBegin(GL_LINES);
					glColor4d(.5f, 0, 1.0f, 1.0f);
					glVertex3d(0,0,0);
					glColor4d(1.0f, 0, 0.5f, .3f);
					glVertex3d((vecLen * system->fluidBox->Vx[system->fluidBox->IX(i, j, k)]),
							   (vecLen * system->fluidBox->Vy[system->fluidBox->IX(i, j, k)]),
							   (vecLen * system->fluidBox->Vz[system->fluidBox->IX(i, j, k)]));
				glEnd();
				//if (system->fluidBox->sphereBnds.find(system->fluidBox->IX(i, j, k)) != system->fluidBox->sphereBnds.end()) {
				//	Eigen::Vector3d normV = system->fluidBox->sphereBnds.find(system->fluidBox->IX(i, j, k))->second->norm;
				//	glPushMatrix();
				//	glColor4f(1.0f, 0.5f, 1.5f, 0.1f);
				//	glutSolidSphere(.25f, 8, 8);
				//	glBegin(GL_LINES);
				//	glColor4d(1.0f, 0,0, .3f);
				//	glVertex3d(0,0,0);
				//	glColor4d(1.0f, 0, 0.5f, 1.0f);
				//	glVertex3d((normV(0)),
				//			   (normV(1)),
				//			   (normV(2)));
				//	glEnd();
				//	glPopMatrix();
				//}
				glPopMatrix();
			}
		}
	}
    glPopMatrix(); 
}//drawFluidVel


void WindowGLDisplay::drawCnstrntLine(vector<std::shared_ptr<myParticle>>& partAra, vector<std::shared_ptr<myConstraint>>& cnstrntAra, vector<vector<int>>& linePointIDS, bool bindCnstr, double cnstrR){
	glColor3f(0.0, 0.0, 0.0);
	int sIdx, eIdx;
	Eigen::Vector3d start, end;
	for each (auto linetup in linePointIDS){
		sIdx = FindPartIDXByID(partAra,linetup[0]);
		//eIdx = (linetup[1] >= 0) ? FindPartIDXByID(partAra,linetup[1]) : MyParticleWorld::systems[MyParticleWorld::curSystem]->getCnstIDXByID(-1 * linetup[1]);//check if 2nd point is actual cnstrt point or center of path constraint
		eIdx = (linetup[1] >= 0) ? FindPartIDXByID(partAra,linetup[1]) : FindCnstrntIDXByID(cnstrntAra,-1 * linetup[1]);//check if 2nd point is actual cnstrt point or center of path constraint

		start = partAra[sIdx]->position[0];
		end = (linetup[1] >= 0) ? partAra[eIdx]->position[0] : cnstrntAra[eIdx]->anchorPoint;//check if 2nd point is actual cnstrt point or center of path constraint
		double lenLine = (start - end).norm();
		if ((!bindCnstr) || ((cnstrR != 0) && (((lenLine/cnstrR) < 1.05 ) || (linetup[1] >= 0)))){
			glBegin(GL_LINES);
				glVertex3d(start[0],start[1], start[2]);
				glVertex3d(end[0], end[1], end[2]);
			glEnd();
		}
	}
}//drawCnstrntLine
void WindowGLDisplay::drawSpring( vector<std::shared_ptr<particleSystem::mySpring>>& springAra){
	vector<int> tmpVec(2, -10);
	//float r;
    float delta_theta = 0.01f;
	Eigen::Vector3d start, end;
	vector<vector<int>> lineCoords;

	glPushMatrix();
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		for(unsigned int sIdx = 0; sIdx < springAra.size(); ++sIdx){
            start = springAra[sIdx]->a->position[0];
			end = springAra[sIdx]->b->position[0];
			//double lenLine = len(start - end);
			glColor3f(0,0,0);
			glBegin(GL_LINES);
				glVertex3d(start[0],start[1], start[2]);
				glVertex3d(end[0], end[1], end[2]);
			glEnd();
		}//for each constraint	
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glPopMatrix();
}//drawCnstrnt
void WindowGLDisplay::drawCnstrnt(vector<std::shared_ptr<particleSystem::myParticle>>& partAra, vector<std::shared_ptr<particleSystem::myConstraint>>& cnstrntAra){
	vector<int> tmpVec(2, -10);
	float r;
    float delta_theta = 0.01f;
	int sIdx, eIdx;
	Eigen::Vector3d start, end;
	vector<vector<int>> lineCoords;

	glPushMatrix();
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glColor3f(0.0, 0.0, 0.0);
		for(unsigned int cidx = 0; cidx < cnstrntAra.size(); ++cidx){
            sIdx = MyParticleWorld::systems[gCurrentScene]->c[cidx]->p1Idx;
            start = partAra[sIdx]->position[0];
			if(cnstrntAra[cidx]->useAnchor){//use anchor means circular path constraint
                    //draw circle
                if(MyParticleWorld::systems[gCurrentScene]->c[cidx]->drawCnstrPath){
				    r = (float)cnstrntAra[cidx]->c_Dist;
				    Eigen::Vector3d loc = cnstrntAra[cidx]->anchorPoint; 
				    glBegin( GL_POLYGON ); // OR GL_LINE_LOOP
				    for( float angle = 0; angle < 2*PI; angle += delta_theta ){
					    glVertex3d( r*cos(angle)+ loc[0], r*sin(angle)+ loc[1], loc[2] );			//draw circle of circular path constrnt
                    }
				    glEnd();
                }

				tmpVec[1] = MyParticleWorld::systems[gCurrentScene]->c[cidx]->ID * -1;
			}//use anchor
			else {//not anchor, draw line regardless of length
				tmpVec[1] = MyParticleWorld::systems[gCurrentScene]->c[cidx]->p2ID;
			}
			eIdx = (tmpVec[1] >= 0) ? FindPartIDXByID(partAra,tmpVec[1]) : FindCnstrntIDXByID(cnstrntAra,-1 * tmpVec[1]);//check if 2nd point is actual cnstrt point or center of path constraint
			end = (tmpVec[1] >= 0) ? partAra[eIdx]->position[0] : cnstrntAra[eIdx]->anchorPoint;//check if 2nd point is actual cnstrt point or center of path constraint
			double cnstrR = (tmpVec[1] >= 0) ? 0 : cnstrntAra[eIdx]->c_Dist;		//radius if constraint is path constraint
			double lenLine = (start - end).norm();
			//if not scene 5 (always draw line)  or current constraint is either closest to a particle (dont' draw multiple path constraint lines to each particle) or not a path constraint
			if ((gCurrentScene != CNSTR_4_JMP) || ((cnstrR != 0) && (((lenLine/cnstrR) < 1.01 ) || (tmpVec[1] >= 0)))){
				//glColor3d(partAra[sIdx]->color[0], partAra[sIdx]->color[1], partAra[sIdx]->color[2]);
				glColor3f(0.0, 0.0, 0.0); 
				glBegin(GL_LINES);
					glVertex3d(start[0],start[1], start[2]);
					glVertex3d(end[0], end[1], end[2]);
				glEnd();
			}
		}//for each constraint	
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glPopMatrix();
}//drawCnstrnt

//draw all particles for a particular scene
void WindowGLDisplay::drawParts(vector<std::shared_ptr<particleSystem::myParticle>>& partAra, double calc_partSize, double d_partSize, bool draw1stPartBlue){
	for(unsigned int i=0; i<partAra.size(); i++){
		if ((draw1stPartBlue) && (i==0)){glColor3f(0.0f, 0.0f, 1.0f);} 
		else {
			glColor3d(partAra[i]->color[0],partAra[i]->color[1],partAra[i]->color[2]);
		}//
		glPushMatrix();
			glTranslated(partAra[i]->position[0][0], partAra[i]->position[0][1], partAra[i]->position[0][2]);
		//	glutSolidSphere(calc_partSize* (log(1 + partAra[i]->mass)), 8, 8);
			glutSolidSphere(calc_partSize, 5, 5);
		glPopMatrix();
	}//for each particle
    //if COM != Eigen::Vector3d(0,0,0), draw it
	if ((MyParticleWorld::systems[gCurrentScene]->partCOM).norm() != 0) {
		glPushMatrix();
    	    glColor3f(0.5f,0.5f,0.5f);
			glTranslated(MyParticleWorld::systems[gCurrentScene]->partCOM[0], MyParticleWorld::systems[gCurrentScene]->partCOM[1], MyParticleWorld::systems[gCurrentScene]->partCOM[2]);
			glutSolidSphere(calc_partSize*1.5, 8, 8);
    	    glColor3f(1.0,1.0,1.0);
			glTranslated(0,0,calc_partSize*1.5);
			glutSolidSphere(calc_partSize, 6, 6);     
        glPopMatrix();

    }
}//drawParts

int WindowGLDisplay::handle(int event)
{  
	Eigen::Vector3d tempPos;
	double xtrans, ytrans;
	std::string lastPress;

	mMouse[0] = fltk::event_x() - (mWinWidth >> 1);
	mMouse[1] = (mWinHeight >> 1) - fltk::event_y();

	int nValidHandle = 0;
	switch(event){

        case fltk::SHORTCUT:{
            int key = fltk::event_key();
            cout<<"Key pressed : ";
            //if(key==32){
            //    cout<<"Space bar"<<endl;
            //    MyParticleWorld::massSprFrcOrVel = !MyParticleWorld::massSprFrcOrVel;
            //    cout<<"Force or Velocity mod on center "<<(MyParticleWorld::massSprFrcOrVel ? "Force" : "Velocity")<<endl;
            //}
            //else 
            {cout<<key<<endl;}

            break;}

		case fltk::PUSH:
			mMouse[0] = fltk::event_x() - (mWinWidth >> 1);
			mMouse[1] = (mWinHeight >> 1) - fltk::event_y();

			MyParticleWorld::flags[MyParticleWorld::msClicked] = true;
			//MyParticleWorld::mouseReleased = false;
			//MyParticleWorld::flags[MyParticleWorld::msDragged] = true;
			MyParticleWorld::msClickVal.pop_back();			
			MyParticleWorld::msClickVal.push_front(mMouse);

			mMouseDown = true;

			mTrackball.StartBall(mMouse[0],mMouse[1]);

			if( fltk::get_key_state(fltk::LeftShiftKey) || fltk::get_key_state(fltk::RightShiftKey) )
				mShiftState = true;
			else
				mShiftState = false;

			if( fltk::get_key_state(fltk::LeftCtrlKey) || fltk::get_key_state(fltk::RightCtrlKey) )
				mControlState = true;
			else
				mControlState = false;

			if( fltk::get_key_state(fltk::LeftAltKey) || fltk::get_key_state(fltk::RightAltKey) )
				mAltState = true;
			else
				mAltState = false;

			if( fltk::get_key_state(fltk::TabKey)  )
				mTabState = true;
			else
				mTabState = false;

			if( fltk::get_key_state('l') || fltk::get_key_state('L') ){
				mLockState = true;
			}else
				mLockState = false;

			if(fltk::event_button() == 1){
				if(fltk::get_key_state(fltk::LeftCtrlKey))
					mWhichClick = 1;
				else if(fltk::get_key_state(fltk::LeftAltKey))
					mWhichClick = 2;
				else 
					mWhichClick = 0;
			}

			if (mWhichClick == 0){	
				if(mShiftState)
					mRotatingXY = true;

				mXstart = mMouse[0];
				mYstart = mMouse[1];

			}else if(fltk::event_button() == 2 || mWhichClick == 1){
				if(mShiftState){
					mTranslating = true;
				}
				else {
					// temp code
					mPicking=true;
					draw();
					mPicking=false;
				}
				mXstart = mMouse[0];
				mYstart = mMouse[1];
			}else if(fltk::event_button() == 3 || mWhichClick == 2){
				if(mShiftState){
					mZooming = true; 
				}
				else {
					// temp code
					mPicking=true;
					draw();
				}
				mXstart = mMouse[0];
				mYstart = mMouse[1]; 

			}
			break;

		case fltk::DRAG:
			if(mRotatingXY){	
				mYangle += (mMouse[0] - mXstart);
				mXangle += (mMouse[1] - mYstart);
				if(mYangle > 360)
					mYangle -= 360;
				else if(mYangle < 0)
					mYangle += 360;
				if(mXangle > 360)
					mXangle -= 360;
				else if(mXangle < 0)
					mXangle +=360;
				mXstart = mMouse[0];
				mYstart = mMouse[1];

				redraw();

				mTrackball.UpdateBall(mMouse[0],mMouse[1]);
			}else if(mZooming){	
				double xVal = mMouse[0] - mXstart;
				double yVal = mMouse[1] - mYstart;
				double shift;

				if(xVal!=0)
					shift = sqrt(xVal * xVal + yVal * yVal) * (xVal / fabs(xVal));
				else
					shift = sqrt(xVal * xVal + yVal * yVal);

				mZoom += shift;
				Eigen::Matrix4d tempRot;
				tempRot.setZero();
				mTrackball.GetCurrentRotation(tempRot);
				tempRot = (tempRot).inverse();
				Eigen::Vector4d tempTrans = tempRot * Eigen::Vector4d(0, 0, shift, 1);
				mFixedTrans[0] += tempTrans[0];
				mFixedTrans[1] += tempTrans[1];
				mFixedTrans[2] += tempTrans[2];

				mXstart = mMouse[0];
				mYstart = mMouse[1];
				redraw();

			}else if(mTranslating){
				mTranslateX += (mMouse[0] - mXstart);
				mTranslateY += (mMouse[1] - mYstart);
				Eigen::Matrix4d tempRot;
				tempRot.setZero();
				mTrackball.GetCurrentRotation(tempRot);
				tempRot = (tempRot).inverse();
				Eigen::Vector4d tempTrans = tempRot * Eigen::Vector4d(mMouse[0] - mXstart, mMouse[1] - mYstart, 0, 1);
				mFixedTrans[0] += tempTrans[0];
				mFixedTrans[1] += tempTrans[1];
				mFixedTrans[2] += tempTrans[2];
				mXstart = mMouse[0];
				mYstart = mMouse[1];
				redraw();

			}else if(mSelectedMark >= 0 && mSelectedMark < nValidHandle){
				xtrans = (mMouse[0] - mXstart)/100.0;
				ytrans = (mMouse[1] - mYstart)/100.0;
				Eigen::Matrix4d tempRot;
				tempRot.setZero();
				Eigen::Matrix4d yRot;
				yRot.setIdentity();
				yRot(0,0) = cos(-M_PI / 2);
				yRot(0,2) = sin(-M_PI / 2);
				yRot(2,0) = -sin(-M_PI / 2);
				yRot(2,2) = cos(-M_PI / 2);

				mTrackball.GetCurrentRotation(tempRot);
				tempRot *= yRot;
				tempRot = (tempRot).inverse();
				Eigen::Vector4d tempTrans = tempRot * Eigen::Vector4d(xtrans, ytrans, 0, 1);
			}
			else if((fltk::event_button() == 2 || mWhichClick == 1)){

			} else {
				MyParticleWorld::flags[MyParticleWorld::msDragged] = true;
			}
			mMouse[0] = fltk::event_x() - (mWinWidth >> 1);
			mMouse[1] = (mWinHeight >> 1) - fltk::event_y();
            mMouse[2] = 0;
			if (MyParticleWorld::flags[MyParticleWorld::msDragged]) {
				MyParticleWorld::msDragVal.pop_back();			
				MyParticleWorld::msDragVal.push_front(mMouse);
			}//most recent value to front
			redraw();
			break;

		case fltk::RELEASE:
			MyParticleWorld::flags[MyParticleWorld::msDragged] = false;
			MyParticleWorld::flags[MyParticleWorld::msClicked] = false;
			MyParticleWorld::flags[MyParticleWorld::msReleased] = true;

			mMouseDown = false;
			mMouse[0] = fltk::event_x() - (mWinWidth >> 1);
			mMouse[1] = (mWinHeight >> 1) - fltk::event_y();
			MyParticleWorld::msRelVal.pop_back();
			MyParticleWorld::msRelVal.push_front(mMouse);

			mPicking = mRotatingXY = mZooming = mTranslating = false;

			if(mBodyTrackBall){
				delete mBodyTrackBall;
				mBodyTrackBall = NULL;
			}

			redraw();
			break;

		default:
			return 0;
	}

	return 1;
}

void WindowGLDisplay::refresh()
{
	redraw();
}

void WindowGLDisplay::resizeWindow(int width, int height)
{
	resize(x(), y(), width, height);
}

void WindowGLDisplay::resize(int x, int y, int w, int h)
{
	mWinWidth = w;
	mWinHeight = h;
	resetCamera();
	GlWindow::resize(x, y, w, h);
}

void WindowGLDisplay::resetCamera()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0,mWinWidth/mWinHeight,.5,600.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable( GL_DEPTH_TEST );
	gluLookAt(mInitEye[0], mInitEye[1], mInitEye[2], mInitLookAt[0], mInitLookAt[1], mInitLookAt[2], mInitUp[0], mInitUp[1], mInitUp[2]);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

Eigen::Vector3d WindowGLDisplay::getNew3DPos(const Eigen::Vector3d& pos) {
	Eigen::Vector3d eye = GetPointEyeToWorld(mInitEye);
	Eigen::Vector3d lookAt = GetPointEyeToWorld(mInitLookAt);

	Eigen::Vector3d temp, newdir;
	GetMouseRay(fltk::event_x(), (mWinHeight-fltk::event_y()), temp, newdir);
	Eigen::Vector3d newpos;
	newpos = eye + newdir*((pos-eye).dot( lookAt-eye))/(newdir.dot(lookAt-eye));
	return newpos;
}
