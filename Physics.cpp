//g++ -o Physics.exe Physics.cpp -L"C:\MinGW\freeglut\include" -lfreeglut -lopengl32 && .\Physics.exe
//.\Physics.exe
#include "GL/glut.h"
#include "GL/freeglut.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <string.h>
#include <initializer_list>
#define PI 3.141592654
#define n 100

class vec{
	public:
	std::vector<double> comps;
	vec(double x, double y, double z): comps({x,y,z}){}
	vec(double x, double y): comps({x,y,0}){}
	vec(){}
	double& operator[](int i){
	return comps[i];
	}
	std::string toString(){
		return std::to_string(comps[0])+","+std::to_string(comps[1])+","+std::to_string(comps[2]);
	}

	vec& operator+=(const vec& r){
		comps[0] += r.comps[0];
		comps[1] += r.comps[1];
		comps[2] += r.comps[2];
		return *this;
	}

	vec& operator-=(const vec& r){
		comps[0] -= r.comps[0];
		comps[1] -= r.comps[1];
		comps[2] -= r.comps[2];
		return *this;
	}
	private:
};


vec operator*(vec v, double d) {return {v.comps[0]*d,v.comps[1]*d,v.comps[2]*d};}

vec operator*(double d, vec v) {return {v.comps[0]*d,v.comps[1]*d,v.comps[2]*d};}

vec operator+(vec l, vec r) {return {l.comps[0]+r.comps[0],l.comps[1]+r.comps[1],l.comps[2]+r.comps[2]};}

vec operator-(vec l, vec r) {return {l.comps[0]-r.comps[0],l.comps[1]-r.comps[1],l.comps[2]-r.comps[2]};}

vec operator-(vec v) {return {-v.comps[0],-v.comps[1],-v.comps[2]};}

vec operator/(vec v, double d) {return {v.comps[0]/d,v.comps[1]/d,v.comps[2]/d};}

double dot(vec v1, vec v2){
	double dot = v1.comps[0]*v2.comps[0]+v1.comps[1]*v2.comps[1]+v1.comps[2]*v2.comps[2];
	return dot;
}

vec normalise(vec v){
	double mag = sqrt(dot(v,v));
	return v/mag;
}

double randvel()
{
    return 6*(rand()/(double(RAND_MAX)+1)-0.5);
}

double dist(vec x2, vec x1){
	double s=sqrt(dot(x2-x1,x2-x1));
	return s;
}

vec Force(vec x2, vec x1,double m1, double m2, double k){
	vec r = x2-x1;
	vec rh= normalise(r);
 	return rh*(k*m1*m2/dot(r,r));
}

vec vnext(vec v,vec x2, vec x1,double m1, double m2, double k,double h){
	vec r = x2-x1;
	vec rh= normalise(r);
	vec k1 = Force(x2,x1,m1,m2,k)/m2;
	vec k2 = Force(x2,x1+h*k1/2,m1,m2,k)/m2;
	vec k3 = Force(x2,x1+h*k2/2,m1,m2,k)/m2;
	vec k4 = Force(x2,x1+h*k3,m1,m2,k)/m2;
	return v+ h*(k1+2.0*k2+2.0*k3+k4)/6.0;
}
vec xnext(vec v,vec x,double h){
	vec k1 = v;
    vec k2 = v + h*k1/2.0;
    vec k3 = v + h*k2/2.0;
    vec k4 = v + h*k3;
    return x + h*(k1+2.0*k2+2.0*k3+k4)/6.0;
}

class Circle{
	public:
    vec x;
	vec v;
    double r;
Circle(vec x,vec v, double r): x(x), r(r), v(v) {}

bool CC_iscol(Circle c2){
	return dot(c2.x-x,c2.x-x)<1.1*(r+c2.r)*(r+c2.r);
}

void resolveCol(Circle &c){
vec d = c.x-x;
vec nor = normalise(d);
x-=nor*(c.r+r-dist(c.x,x))/2;
c.x+= nor*(c.r+r-dist(c.x,x))/2;
//nor = normalise(d);
vec relv=c.v-v;
double sc = relv[0]*nor[0]+relv[1]*nor[1];
if(sc>0) return;
v+=sc*nor;
c.v-= sc*nor;
}
void simulate(vec acc){
	v += acc*(1.0/100.0);
	x += v;
}

};

class GBall{
		public:
    vec x;
	vec v;
    double r;
	double m;
GBall(vec x,vec v, double r,double m): x(x), r(r), v(v), m(m)
{}
bool GB_iscol(GBall c2){
	return dot(c2.x-x,c2.x-x)<1.1*(r+c2.r)*(r+c2.r);
}

void resolveCol(GBall &c){
vec d = c.x-x;
vec nor = normalise(d);
x-=nor*(c.r+r-dist(c.x,x))/2;
c.x+= nor*(c.r+r-dist(c.x,x))/2;
//nor = normalise(d);
vec relv=c.v-v;
double sc = relv[0]*nor[0]+relv[1]*nor[1];
if(sc>0) return;
v+=sc*nor;
c.v-= sc*nor;
}
void resolveInter(GBall &c){
	v = vnext(v,c.x,x,c.m,m,10,0.001);
	c.v = vnext(c.v,x,c.x,m,c.m,10,0.001);
}
void simulate(){
	x = xnext(v,x,0.001);
	//x+=v;
}

};

struct Line{
    vec xi;
    vec xf;
Line(vec xi,vec xf): xi(xi), xf(xf) {}
};

std::vector<Line> lines = {};
std::vector<Circle> circles ={};
std::vector<GBall> Balls ={};

void addcircle(double x1,double x2,double v1,double v2,double r){
	Circle c ={{x1,x2},{v1,v2},r};
	circles.push_back(c);
}

void addball(double x1,double x2,double v1,double v2,double r,double m){
	GBall c ={{x1,x2},{v1,v2},r,m};
	Balls.push_back(c);
}

void randpart(double x1,double x2,double r){
	Circle c = {{x1,x2},{randvel(),randvel()},r};
	circles.push_back(c);
}

void linegen(double x1,double x2,double x3,double x4){
	Line l = {{x1,x2},{x3,x4}};
	lines.push_back(l);
}

void init_gas(int m){
	const double r=4,ip=-225,fp=225;
	double ds = (fp-ip)/double(m);
	for(int i=0; i<m;i++){
		for(int j=0;j<m;j++){
			randpart(ip+i*ds,ip+j*ds,r);
		}
	}
}

double linedist(vec vi, vec vf, vec x){
    vec v = vf-vi;
    vec vp= {-v[1],v[0]};
    double nu =std::abs(dot(vp,(x-vi)));
    double d = sqrt(dot(vp,vp));
    return nu/d;
}

bool LC_iscol(Circle c, Line l){
	bool first = linedist(l.xi,l.xf,c.x)<c.r;
	return first;
}

double depthPLC(Circle c, Line l){
	return c.r-linedist(l.xi,l.xf,c.x);
}


vec dpLC(Circle c, Line l){
	vec v = l.xf-l.xi;
	vec vp = {-v[1],v[0]};
    vec nor = normalise(vp);
	double sc = -2.*dot(c.v,nor);
	return sc*nor;
}

vec correctLC(Circle c, Line l){
vec v = l.xf-l.xi;
vec vp = {-v[1],v[0]};
vec nor = normalise(vp);
return nor*depthPLC(c,l);
}

void changeViewPort(int w, int h)
{
	glViewport(0, 0, w, h);
}


void circle(double x, double y, double r){
    glBegin(GL_LINE_STRIP);
	for (int j=0;j<n;j++){
		glColor3f(0.0f, 0.0f, 0.0f);
	glVertex2f(x+r*cos(j*2*PI/n),y+r*sin(j*2*PI/n));
	}
	glEnd();
}
void circler(double x, double y, double r){
    glBegin(GL_LINE_STRIP);
	for (int j=0;j<n;j++){
		glColor3f(0.0f, 1.0f, 0.0f);
	glVertex2f(x+r*cos(j*2*PI/n),y+r*sin(j*2*PI/n));
	}
	glEnd();
}

void line(float xc1,float yc1,float xc2,float yc2)
{
	glBegin(GL_LINES);
	glColor3f(0.4f, 0.4f, 0.4f);
	glVertex2f(xc1,yc1);
	glVertex2f(xc2,yc2);
	glEnd();
}
void linew(float xc1,float yc1,float xc2,float yc2)
{
	glBegin(GL_LINES);
	glColor3f(0.0f, 0.0f, 0.0f);
	glVertex2f(xc1,yc1);
	glVertex2f(xc2,yc2);
	glEnd();
}

void lineb(float xc1,float yc1,float xc2,float yc2)
{
	glBegin(GL_LINES);
	glColor3f(0.0f, 0.0f, 0.6f);
	glVertex2f(xc1,yc1);
	glVertex2f(xc2,yc2);
	glEnd();
}

void grid (float l,int num){
	for(int i=0;i<num;i++)
{
line(i*l-300,300,i*l-300,-300);
}
	for(int i=0;i<num;i++)
{
line(300,i*l-300,-300,i*l-300);
}

}

void rendervels(){
	const double r=2;
	double maxv =0;
for(int i=0;i<circles.size();i++){
		Circle &c = circles[i];
		double spdi = sqrt(dot(c.v,c.v));
		if (spdi>=maxv){
			maxv=spdi;
		}
	}	
	double sc= 500.0/maxv;
for(int i=0;i<circles.size();i++){
		Circle &c = circles[i];
		double speed = sqrt(dot(c.v,c.v));
		circle(sc*speed-250.,-275.,r);
	}
}

void rendercm(){
	double s=0;
	vec cm ={0,0};
	for(int i=0;i<Balls.size();i++){
		GBall &c = Balls[i];
		cm=cm+c.x*c.m;
		s+=c.m;
	}
	cm= cm/s;
	circler(cm[0],cm[1],3);
}

void render()
{
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH, GL_NICEST);
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH, GL_NICEST);
	glClearColor(1.0f,1.0f,1.0f,1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0,0.0,0.0); 
	glOrtho(-400.0,400.0, -320.0,320.0, 0, 1);
	//desplaza todo
	//glPushMatrix(), glPopMatrix() Transforma el espacio relativamente a la figura
	//grid(15,40);
    //grid(10,60);
    //linew(0.,-300.,0.,300.);
    //linew(-300.,-250.,300.,-250.);

	for(int i=0;i<circles.size();i++){
		Circle &c = circles[i];
		circle(c.x[0],c.x[1],c.r);
	}
	for(int i=0;i<Balls.size();i++){
		GBall &c = Balls[i];
		circle(c.x[0],c.x[1],c.r);
	}

	for(int i=0;i<lines.size();i++){
		Line &l = lines[i];
		lineb(l.xi[0], l.xi[1],l.xf[0], l.xf[1]);
	}
	//rendervels();
	//rendercm();
	glutSwapBuffers();
}

void SimCircles(vec g){
for(int i=0;i<circles.size();i++){
		Circle &c = circles[i];
		c.simulate(g);
		for(int j=i+1;j<circles.size();j++){
		Circle &c2 = circles[j];
		if(c.CC_iscol(c2)){
			c.resolveCol(c2);
		}
		}
		for(int i=0;i<lines.size();i++){
		Line &l = lines[i];
		if (LC_iscol(c,l)){
			c.x+=correctLC(c,l);
			c.v+=dpLC(c,l);	
			}
		}
	}
}

void SimForceB(){
	for(int i=0;i<Balls.size();i++){
		GBall &c = Balls[i];
		c.simulate();
		for(int j=i+1;j<Balls.size();j++){
		GBall &c2 = Balls[j];
		if(!c.GB_iscol(c2)){
			c.resolveInter(c2);
		}
		else c.resolveCol(c2);
		}
		
	}
}

void ballFromString(std::string str){

}

void time (int v){
	//vec g = {0.,-0.01};
	//SimCircles(g);
	SimForceB();
	render();
	glutPostRedisplay();
	glutTimerFunc(1,time,0);
}

int main(int argc, char* argv[]) {
	linegen(-380,-300,380,-300);
	linegen(-380,300,380,300);
	linegen(-380,-300,-380,300);
	linegen(380,300,380,-300);
	//init_gas(15);
	addball(0,-150,0,1,5,1000);
	addball(50,-150,0,15,5,1);

	
	// Initialize GLUT
	glutInit(&argc, argv);

	// Set up some memory buffers for our display
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	// Set the window size
	glutInitWindowSize(800, 640);
	// Create the window with the title "Hello,GL"
	glutCreateWindow("GLarmad");
	// Bind the two functions (above) to respond when necessary
	glutReshapeFunc(changeViewPort);
	glutDisplayFunc(render);
	// Very important!  This initializes the entry points in the OpenGL driver so we can 
	// call all the functions in the API.
	glutTimerFunc(1,time,0);
	glutMainLoop();
	
	return 0;
}