#ifdef WIN32
#include <windows.h>
#endif
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>

#include <GL/gl.h>
#include <GL/glu.h>
#include "../../external/glut/GL/glut.h"
#include <fltk/visual.h>

extern int gCurrentScene;
extern bool gReset;

#include "myGlobConsts.h"
#include "myParticle.h"
#include "myConstraint.h"
#include "mySystem.h"

#include "WindowGLDisplay.h"
#include "ControlPanel.h"
#include "MyParticleWorld.h"
#include "ParticleSystemUI.h"

ParticleSystemUI* mUI;

using namespace std;
using namespace particleSystem;

void onTimer(void *){
	double howFast = .001;
	if (!MyParticleWorld::systems[gCurrentScene]->flags[mySystem::pauseSim]) {
		if( mUI->mControl->isPlaying()){//playback frames
            recordFrames();
			mUI->mControl->advanceFrame();
		}
		if (gReset){
			gReset = false;
			mUI->Show();
		}
		if(mUI->mSim_but->value()){	
			MyParticleWorld::handleTimeStep(gCurrentScene);
		}
		howFast = mUI->mControl->getSpeed();
    } else {//sim is paused
		MyParticleWorld::handlePause(gCurrentScene);
	}
	mUI->mGLWindow->refresh();
	fltk::add_timeout( howFast, onTimer );	
}

int main(int argc, char ** argv){

	// inits
	gReset = false;
	// ---------------------------------------------
	mUI = new ParticleSystemUI;
	// ---------------------------------------------
	for (int i = 0; i < numCurrSystems; ++i) {
		MyParticleWorld::initScene(i, globDeltaT, RK4_G, 0, 0, 0);			
	}
	gCurrentScene = MSPR_MTN_PROJ;		//start with 30 sided die

    // show UI
	fltk::visual(fltk::DOUBLE_BUFFER|fltk::INDEXED_COLOR);
	mUI->Show();
	fltk::add_timeout( globDeltaT, onTimer );
	return fltk::run();
}
