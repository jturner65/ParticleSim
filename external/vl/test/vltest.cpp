#include "vl\VLfd.h"

void Test(Vec2f a, Vec3f b, Vec4f c)
{
	cout << a << b << c << endl;
}

void AssignmentTest()
{
	Vec2f v2 = vl_one;
	Vec3f v3 = vl_one;
	Vec4f v4 = vl_one;
	Vecf  vv(10, vl_one);
	Mat2d m2 = vl_I;
	Mat3d m3 = vl_I;
	Mat4d m4 = vl_one;
	Matd  mm(10,6); mm = vl_I;

	cout << v2 << v3 << v4 << vv << endl;
	cout << m2 << m3 << m4 << mm << endl;
		
	Test(vl_zero, vl_y, vl_one);
}

void main()
{
	Mat4d M(1.0, 2.0, 3.0, 4.0, 
		5.0, 6.0, 7.0, 8.0, 
		9.0, 10.0, 11.0, 12.0, 
		13.0, 14.0, 15.0, 16.0);
	Vec3f v(1.0, 2.0, 3.0);

	// homogeneous multiply
	
	cout << proj(M * Vec4f(v, 1.0)) << endl;
	cout << xform(M, v) << endl;
	
	AssignmentTest();
}
