// generated by Fast Light User Interface Designer (fluid) version 2.0100

#ifndef mandelbrot_ui_h
#define mandelbrot_ui_h
#include "mandelbrot.h"
#include <fltk/Window.h>
#include <fltk/FloatInput.h>
#include <fltk/Slider.h>
#include <fltk/InvisibleBox.h>

class Drawing_Window  {
public:
  void make_window();
  fltk::Window *window;
    Drawing_Area *d;
    fltk::FloatInput *x_input;
private:
    inline void cb_x_input_i(fltk::FloatInput*, void*);
    static void cb_x_input(fltk::FloatInput*, void*);
public:
    fltk::FloatInput *y_input;
private:
    inline void cb_y_input_i(fltk::FloatInput*, void*);
    static void cb_y_input(fltk::FloatInput*, void*);
public:
    fltk::FloatInput *w_input;
private:
    inline void cb_w_input_i(fltk::FloatInput*, void*);
    static void cb_w_input(fltk::FloatInput*, void*);
    inline void cb_brightness_i(fltk::Slider*, void*);
    static void cb_brightness(fltk::Slider*, void*);
    inline void cb_iterations_i(fltk::Slider*, void*);
    static void cb_iterations(fltk::Slider*, void*);
public:
  void update_label();
};
#endif
