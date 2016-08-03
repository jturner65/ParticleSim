//
// "$Id: setcolor.cxx 5560 2006-12-20 03:29:32Z spitzak $"
//
// MacOS color functions for the Fast Light Tool Kit (FLTK).
//
// Copyright 1998-2006 by Bill Spitzak and others.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA.
//
// Please report all bugs and problems to "fltk-bugs@fltk.org".
//

// we need to remember some settings for the current context

namespace fltk {

static float quartz_line_width_ = 1.0f;
static enum CGLineCap quartz_line_cap_ = kCGLineCapButt;
static enum CGLineJoin quartz_line_join_ = kCGLineJoinMiter;
static float *quartz_line_pattern = 0;
static int quartz_line_pattern_size = 0;

void restore_quartz_line_style() {
  CGContextSetLineWidth(quartz_gc, quartz_line_width_);
  CGContextSetLineCap(quartz_gc, quartz_line_cap_);
  CGContextSetLineJoin(quartz_gc, quartz_line_join_);
  CGContextSetLineDash(quartz_gc, 0,
      quartz_line_pattern, quartz_line_pattern_size);
}

} // end of namespace

// Because carbon has a 'current color' in the drawing context this
// is really simple.

void fltk::setcolor(Color i) {
  current_color_ = i;
  // get fltk indexed color:
  if (!(i & 0xFFFFFF00)) i = cmap[i];
  // get the individual colors and put into Mac color structure:
  if (!quartz_gc) return; // no context yet? We will assign the color later.
  uchar r, g, b;
  r = i>>24;
  g = i>>16;
  b = i>> 8;
#if USE_CAIRO
  cairo_set_source_rgb(cc, r/255.0, g/255.0, b/255.0);
#endif
  float fr = r/255.0f;
  float fg = g/255.0f;
  float fb = b/255.0f;
  CGContextSetRGBFillColor(quartz_gc, fr, fg, fb, 1.0f);
  CGContextSetRGBStrokeColor(quartz_gc, fr, fg, fb, 1.0f);
}

// Used by setcolor_index
static void free_color(Color) {}

static enum CGLineCap Cap[4] = {
  kCGLineCapButt, kCGLineCapButt, kCGLineCapRound, kCGLineCapSquare
};

static enum CGLineJoin Join[4] = {
  kCGLineJoinMiter, kCGLineJoinMiter, kCGLineJoinRound, kCGLineJoinBevel
};

void fltk::line_style(int style, float  width, const char* dashes) {
  line_style_ = style;
  line_width_ = width;
  line_dashes_ = dashes;
  quartz_line_width_ = width ? width : 1;
  quartz_line_cap_ = Cap[(style>>8)&3];
  quartz_line_join_ = Join[(style>>12)&3];
  const char *d = dashes;
  static float pattern[16];
  if (d && *d) {
    float *p = pattern;
    while (*d) { *p++ = (float)*d++;}
    quartz_line_pattern = pattern;
    quartz_line_pattern_size = d-dashes;
  } else if (style & 0xff) {
    float dash, dot, gap;
    // adjust lengths to account for cap:
    if (style & 0x200) {
      dash = 2*quartz_line_width_;
      dot = 0;
      gap = 2*quartz_line_width_;
    } else {
      dash = 3*quartz_line_width_;
      dot = gap = quartz_line_width_;
    }
    float *p = pattern;
    switch (style & 0xff) {
    case DASH:       *p++ = dash; *p++ = gap; break;
    case DOT:        *p++ = dot; *p++ = gap; break;
    case DASHDOT:    *p++ = dash; *p++ = gap; *p++ = dot; *p++ = gap; break;
    case DASHDOTDOT: *p++ = dash; *p++ = gap; *p++ = dot; *p++ = gap;
                     *p++ = dot; *p++ = gap; break;
    }
    quartz_line_pattern_size = p-pattern;
    quartz_line_pattern = pattern;
  } else {
    quartz_line_pattern = 0;
    quartz_line_pattern_size = 0;
  }
  restore_quartz_line_style();
#if USE_CAIRO
  cairo_set_line_width(cc, width ? width : 1.0);
  int c = (style>>8)&3; if (c) c--;
  cairo_set_line_cap(cc, (cairo_line_cap_t)c);
  int j = (style>>12)&3; if (j) j--;
  cairo_set_line_join(cc, (cairo_line_join_t)j);
  int ndashes = quartz_line_pattern_size;
  if (ndashes) {
    double dash[20];
    for (int i = 0; i < ndashes; i++) dash[i] = (double) dashes[i];
    cairo_set_dash(cc, dash, ndashes, 0);
  } else {
    cairo_set_dash(cc, 0, 0, 0);
  }
#endif
}

//
// End of "$Id: setcolor.cxx 5560 2006-12-20 03:29:32Z spitzak $".
//
