#include "glFuncs.h"
#include <iostream>
#if USEFLTK
#include <fltk/file_chooser.h>
#include <fltk/Style.h>
#endif

using namespace std;

bool ScreenShot(int w, int h, char *fname, bool _antialias) {
	// make sure OpenGL context is current
	//make_current();

	// read the pixels
	int numPixels = w*h;
	unsigned char *pixels = new unsigned char[numPixels*3*sizeof(unsigned char)];
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	if(_antialias) glReadBuffer(GL_ACCUM);
	else glReadBuffer(GL_FRONT);
	glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	// swap red and blue, because TGA format is stupid
	int i;
	for (i=0; i < numPixels; i++) {
		pixels[i * 3 + 0] ^= pixels[i * 3 + 2];
		pixels[i * 3 + 2] ^= pixels[i * 3 + 0];
		pixels[i * 3 + 0] ^= pixels[i * 3 + 2];
	}

	// get file name
#if USEFLTK
	if (fname == NULL) {
		fname = (char *)fltk::file_chooser("L", "*.tga", NULL);
	}
#endif
	if (fname == NULL)
		return false;

	// open the file
	FILE *fptr;
	fptr = fopen(fname, "wb");
	if (fptr == NULL) {
		//    PrintOnScreen("Failed to open this file");
		return false;
	}

	// create tga header
	putc(0,fptr);
	putc(0,fptr);
	putc(2,fptr);                         // uncompressed RGB
	putc(0,fptr); putc(0,fptr);
	putc(0,fptr); putc(0,fptr);
	putc(0,fptr);
	putc(0,fptr); putc(0,fptr);           // X origin
	putc(0,fptr); putc(0,fptr);           // y origin
	putc((w & 0x00FF),fptr);
	putc((w & 0xFF00) / 256,fptr);
	putc((h & 0x00FF),fptr);
	putc((h & 0xFF00) / 256,fptr);
	putc(24,fptr);                        // 24 bit bitmap
	putc(0,fptr);

	// write the data
	fwrite(pixels, w*h*3*sizeof(char), 1, fptr);
	fclose(fptr);

	delete []pixels;

	cout << fname << " generated" << endl;
	return true;
}

void DrawStringOnScreen2(double x, double y, char *s)
{ // draws text on the screen
	GLint oldMode;
	glGetIntegerv(GL_MATRIX_MODE, &oldMode);
	glMatrixMode( GL_PROJECTION );

	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D( 0.0, 1.0, 0.0, 1.0 );

	glMatrixMode( GL_MODELVIEW );
	glPushMatrix();
	glLoadIdentity();
	fltk::glsetfont(fltk::HELVETICA, 12);
	fltk::gldrawtext(s, (float)x, (float)y);
	glPopMatrix();

	glMatrixMode( GL_PROJECTION );
	glPopMatrix();
	glMatrixMode(oldMode);
}

void DrawStringOnScreen(float x, float y, void *font, std::string s){
	glRasterPos2f(x, y);                             // set position to start drawing fonts
	for (unsigned int c=0; c < s.length(); c++)
		glutBitmapCharacter(font, s.at(c) );  // draw the character to the screen
}

void DrawArrow(const Eigen::Vector3d& pt, const Eigen::Vector3d& dir, double length, double thickness, double arrowthickness) {
	Eigen::Vector3d normDir = dir.normalized();

	if(arrowthickness==-1) arrowthickness=2*thickness;
	double arrowlength = 2*arrowthickness;

	GLUquadricObj *c;
	c = gluNewQuadric();
	gluQuadricDrawStyle(c, GLU_FILL);
	gluQuadricNormals(c, GLU_SMOOTH);

	glPushMatrix();
	glTranslated(pt[0], pt[1], pt[2]);
	glRotated(acos(normDir[2])*180/M_PI, -normDir[1], normDir[0], 0);
	gluCylinder(c, thickness, thickness, length-arrowlength, 16, 16);

	// arrowhed
	glPushMatrix();
	glTranslated(0, 0, length-arrowlength);
	gluCylinder(c, arrowthickness, 0.0, arrowlength, 10, 10);
	glPopMatrix();

	glPopMatrix();

	gluDeleteQuadric(c);
}

void GetMouseRay(double x, double y, Eigen::Vector3d& pos, Eigen::Vector3d& dir) {
	GLdouble px,py,pz;
 	GLdouble modelmat[16];
	GLdouble projmat[16];
	GLint viewport[4];

	glGetDoublev( GL_MODELVIEW_MATRIX, modelmat );
	glGetDoublev( GL_PROJECTION_MATRIX, projmat );
	glGetIntegerv( GL_VIEWPORT, viewport );

	gluUnProject(x, y, 0.0, modelmat, projmat, viewport, &px, &py, &pz);

	pos[0] = px;
	pos[1] = py;
	pos[2] = pz;

	gluUnProject(x, y, 1.0, modelmat, projmat, viewport, &px, &py, &pz);

	dir[0] = px - pos[0];
	dir[1] = py - pos[1];
	dir[2] = pz - pos[2];

	dir.normalize();
}

Eigen::Vector3d GetCameraCenterWorld() {
	GLdouble projmat[16];
	GLdouble modelmat[16];

	glGetDoublev( GL_PROJECTION_MATRIX, projmat);
	glGetDoublev( GL_MODELVIEW_MATRIX, modelmat);
	Eigen::Matrix4d InvModelMat;
	Eigen::Matrix4d InvProjMat;
	InvModelMat.setZero();
	InvProjMat.setZero();
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			InvModelMat(i,j)=modelmat[4*i+j];
			InvProjMat(i,j)=projmat[4*i+j];
			//printf("%lf ", InvModelMat[i][j]);
		}
		//printf("\n");
	}
	InvModelMat=(InvModelMat*InvProjMat).transpose();
	InvModelMat = (InvModelMat.inverse());
	Eigen::Vector4d cc(0, 0, 0, 1);
	cc=InvModelMat*cc;
	if(cc.norm()!=0) cc=cc/cc[3];
	return Eigen::Vector3d(cc[0], cc[1], cc[2]);
}

Eigen::Vector3d GetPointEyeToWorld(const Eigen::Vector3d& _pt) {
	GLdouble modelmat[16];

	glGetDoublev( GL_MODELVIEW_MATRIX, modelmat);
	Eigen::Matrix4d InvModelMat;
	InvModelMat.setZero();
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			InvModelMat(i,j)=modelmat[4*i+j];
		}
	}
	InvModelMat=(InvModelMat.transpose());
	InvModelMat=(InvModelMat.inverse());
	Eigen::Vector4d cc( _pt(0),_pt(1), _pt(2),1);
	cc=InvModelMat*cc;
	if(cc.norm()!=0) cc=cc/cc[3];
	return Eigen::Vector3d(cc[0], cc[1], cc[2]);
}


void accFrustum(GLdouble _left, GLdouble _right, GLdouble _bottom,
				GLdouble _top, GLdouble _near, GLdouble _far, GLdouble _pixdx, 
				GLdouble _pixdy, GLdouble _eyedx, GLdouble _eyedy, 
				GLdouble _focus)
{
	GLdouble xwsize, ywsize; 
	GLdouble dx, dy;
	GLint viewport[4];

	glGetIntegerv (GL_VIEWPORT, viewport);

	xwsize = _right - _left;
	ywsize = _top - _bottom;
	dx = -(_pixdx*xwsize/(GLdouble)viewport[2] + _eyedx*_near/_focus);
	dy = -(_pixdy*ywsize/(GLdouble)viewport[3] + _eyedy*_near/_focus);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum (_left + dx, _right + dx, _bottom + dy, _top + dy, 
		_near, _far);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslated (-_eyedx, -_eyedy, 0.0);
	glMatrixMode(GL_PROJECTION);
}

void accPerspective(GLdouble _fovy, GLdouble _aspect, 
					GLdouble _near, GLdouble _far, GLdouble _pixdx, GLdouble _pixdy, 
					GLdouble _eyedx, GLdouble _eyedy, GLdouble _focus)
{
	GLdouble fov2,left,right,bottom,top;
	fov2 = ((_fovy*M_PI) / 180.0) / 2.0;

	top = _near / (cosf(fov2) / sinf(fov2));
	bottom = -top;
	right = top * _aspect;
	left = -right;

	accFrustum (left, right, bottom, top, _near, _far,
		_pixdx, _pixdy, _eyedx, _eyedy, _focus);
}


int LoadBMP(char *filename, Image *image)
{
	FILE *file;
	unsigned long size;        // size of the image in bytes.
	unsigned long i;           // standard counter.
	unsigned short int planes; // number of planes in image (must be 1) 
	unsigned short int bpp;    // number of bits per pixel (must be 24)
	char temp;                 // temporary color storage for bgr-rgb conversion.

	// The job of this code is to
	// load in a bitmap file. If the file doesn't exist NULL is sent back
	// meaning the texture couldn't be loaded. Before I start explaining the
	// code there are a few VERY important things you need to know about the
	// images you plan to use as textures. The image height and width MUST be a
	// power of 2. The width and height must be at least 64 pixels, and for
	// compatibility reasons, shouldn't be more than 256 pixels. If the image
	// you want to use is not 64, 128 or 256 pixels on the width or height,
	// resize it in an art program. There are ways around this limitation, but
	// for now we'll just stick to standard texture sizes. 

	// make sure the file is there.

	if ((file = fopen(filename, "rb"))==NULL) {
		printf("File Not Found : %s\n",filename);
		return 0;
	}

	// seek through the bmp header, up to the width/height:
	fseek(file, 18, SEEK_CUR);

	// No 100% errorchecking anymore!!!

	// read the width
	fread(&image->sizeX, sizeof(int), 1, file);
	printf("Width of %s: %lu\n", filename, image->sizeX);

	// read the height
	fread(&image->sizeY, sizeof(int), 1, file);
	printf("Height of %s: %lu\n", filename, image->sizeY);

	// calculate the size (assuming 24 bits or 3 bytes per pixel).
	size = image->sizeX * image->sizeY * 4;

	// read the planes
	fread(&planes, sizeof(short), 1, file);
	if (planes != 1) {
		printf("Planes from %s is not 1: %u\n", filename, planes);
		return 0;
	}

	// read the bpp
	fread(&bpp, sizeof(short), 1, file);
	if (bpp != 24) {
		printf("Bpp from %s is not 24: %u\n", filename, bpp);
		return 0;
	}

	// seek past the rest of the bitmap header.
	fseek(file, 24, SEEK_CUR);

	// read the data. 
	image->pixels = (unsigned char *) malloc(size);
	if (image->pixels == NULL) {
		printf("Error allocating memory for color-corrected image data");
		return 0;
	}

	for (i=0;i<size;i+=4) { // reverse all of the colors. (bgr -> rgb)
		if (fread(&image->pixels[i], 3, 1, file) != 1) {
			printf("Error reading image data from %s.\n", filename);
			return 0;
		}
		temp = image->pixels[i];
		image->pixels[i] = image->pixels[i+2];
		image->pixels[i+2] = temp;
		image->pixels[i+3] = 255;
	}

	// we're done.
	return 1;
}

// Load Bitmaps And Convert To Textures
unsigned int load_texture(char *filename)
{
	// Load Texture
	Image image_data;
	unsigned int   texture_num;

	// allocate space for texture
	if (!LoadBMP(filename, &image_data)) {
		std::exit(1);
	}

	// Create Texture   
	glGenTextures(1, &texture_num);
	// 2d texture (x and y size)
	glBindTexture(GL_TEXTURE_2D, texture_num);

	// scale linearly when image bigger than texture
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR); 

	// scale linearly when image smalled than texture
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR); 

	// 2d texture, level of detail 0 (normal), 3 components 
	// (red, green, blue), x size from image, y size from image, 
	// border 0 (normal), rgb color data, unsigned byte data, and 
	// finally the data itself.

	glTexImage2D(GL_TEXTURE_2D, 0, 3, image_data.sizeX, image_data.sizeY, 0, 
		GL_RGBA, GL_UNSIGNED_BYTE, image_data.pixels);

	free(image_data.pixels);
	return texture_num;
}
