#ifndef __UICALLBACK_H__
#define __UICALLBACK_H__

#include <fltk/run.h>
#include <fltk/Widget.h>


class fltk::Widget;

void Exit_cb(fltk::Widget *o, void *v);
void StillShot_cb(fltk::Widget *o, void *v);
void Sim_cb(fltk::Widget *o, void *v);
void Pause_cb(fltk::Widget *o, void *v);
void Switch_cb(fltk::Widget *o, void *v);

void recordFrames();

#endif
