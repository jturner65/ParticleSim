// generated by Fast Light User Interface Designer (fluid) version 2.1000

#ifndef ControlPanel_h
#define ControlPanel_h
#include <fltk/DoubleBufferWindow.h>
#include <fltk/ValueSlider.h>
#include <fltk/Button.h>
#include <fltk/ValueInput.h>

class ControlPanel  {
public:
  void setRange(int _maxFrame);
  ControlPanel();
  fltk::DoubleBufferWindow *mControlWindow;
    fltk::ValueSlider *mFrame_sli;
private:
    inline void cb_mFrame_sli_i(fltk::ValueSlider*, void*);
    static void cb_mFrame_sli(fltk::ValueSlider*, void*);
public:
    fltk::ValueSlider *mBegin_sli;
private:
    inline void cb_mBegin_sli_i(fltk::ValueSlider*, void*);
    static void cb_mBegin_sli(fltk::ValueSlider*, void*);
public:
    fltk::ValueSlider *mEnd_sli;
private:
    inline void cb_mEnd_sli_i(fltk::ValueSlider*, void*);
    static void cb_mEnd_sli(fltk::ValueSlider*, void*);
public:
    fltk::Button *mLoop_but;
    fltk::Button *mPlay_but;
private:
    inline void cb_mPlay_but_i(fltk::Button*, void*);
    static void cb_mPlay_but(fltk::Button*, void*);
public:
    fltk::ValueInput *mFrameStep_inp;
    fltk::Button *mForward_but;
private:
    inline void cb_mForward_but_i(fltk::Button*, void*);
    static void cb_mForward_but(fltk::Button*, void*);
public:
    fltk::Button *mBackward_but;
private:
    inline void cb_mBackward_but_i(fltk::Button*, void*);
    static void cb_mBackward_but(fltk::Button*, void*);

public://mPartCount_inp
    fltk::ValueInput *mPartCount_inp;
	int getCurPartCount();
    fltk::Button *mForward_butP;
private:
    inline void cb_mForward_butP_i(fltk::Button*, void*);
    static void cb_mForward_butP(fltk::Button*, void*);
public:
    fltk::Button *mBackward_butP;
private:
    inline void cb_mBackward_butP_i(fltk::Button*, void*);
    static void cb_mBackward_butP(fltk::Button*, void*);


public:
    fltk::ValueSlider *mSpeed_rol;
  void advanceFrame(int _step=1);
  void decrementFrame(int _step=1);
  int getCurrentFrame();
  void setCurrentFrame(int _f);
  int getMaxFrame();
  int getBegin();
  int getEnd();
  double getSpeed();
  bool isLooping();
  bool isPlaying();
private:
  bool mPlaying;
};
#endif
