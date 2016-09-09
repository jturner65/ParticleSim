#include "main.h"
#include "UICallback.h"
#include "ControlPanel.h"
#include "WindowGLDisplay.h"
#include "glFuncs.h"
#include "MyParticleWorld.h"
#include "mySystem.h"
#include "ParticleSystemUI.h"

#include <vector>
#include <fltk/file_chooser.h>

using namespace std;

bool gReset;
int gCurrentScene;

extern ParticleSystemUI* mUI;

void Exit_cb(fltk::Widget *o, void *v)
{
	exit(0);
}

void StillShot_cb(fltk::Widget *o, void *v)
{
	// take a still shot
	static int count=0;
	char fileBase[32]="ScreenShot";
	char fileName[64];
	sprintf(fileName, "%s%.4d.tga", fileBase, count++); 
	ScreenShot(mUI->mGLWindow->w(), mUI->mGLWindow->h(), fileName);
}

void recordFrames(){
	// saving the dofs
	if(mUI->mSim_but->value()){

	}

	// saving the images
	if(mUI->mRecordMotion_but->value() && (mUI->mSim_but->value() || mUI->mControl->isPlaying())){
		static int count=0;
		char fileBase[32]="Capture";
		char fileName[64];
		sprintf(fileName, "%s%.4d.tga", fileBase, count++); 
		ScreenShot(mUI->mGLWindow->w(), mUI->mGLWindow->h(), fileName, mUI->mGLWindow->mAntiAlias);
	}
}

void Sim_cb(fltk::Widget *o, void *v){
	if(mUI->mSim_but->value()){
		mUI->mControlGrp->deactivate();
	}
	else {
		mUI->mControlGrp->activate();
		int numFrames = MyParticleWorld::systemParticleHolder.size();
		mUI->mControl->setRange(numFrames-1);
	}
}

void Pause_cb(fltk::Widget *o, void *v){
	if(mUI->mPause_but->value()){
		MyParticleWorld::systems[gCurrentScene]->flags[mySystem::pauseSim] = true;
	}
	else {
		MyParticleWorld::systems[gCurrentScene]->flags[mySystem::pauseSim] = false;
	}
}

void Switch_cb(fltk::Widget *o, void *v){
	gCurrentScene = (gCurrentScene + 1) % numCurrSystems;
	//gCurrentScene++;
	//if (gCurrentScene >= numCurrSystems){gCurrentScene = 0;}
	gReset = true;
    mUI->mControl->setRange(0);
}
