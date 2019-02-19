#include <iostream>
#include <QMessageBox>

#include "gl_view.h"
#include "gl_utils.h"
#include "height_calculator.h"
#include "main_dialog.h"
#include "world.h"
#include "vec3.h"

using namespace std;
using namespace Eigen;
using namespace farn;

const double rho = 5.96e3;

GlView::GlView(MainDialog *parent):
  QGLWidget(parent),
  pHeightCalculator(new HeightCalculator),
  shiftX(0),
  shiftY(0),
  pMainDialog(parent),
  pGraph(new GraphWidget())
  {
    xang = 0;
    yang = 0;
    zang = 0;
  }


void
GlView::initializeGL(){
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glEnable(GL_BLEND);

  glEnable(GL_SMOOTH);
  //glEnable(GL_MULTISAMPLE);

  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

  glEnable(GL_COLOR_MATERIAL);

  glEnable(GL_LIGHTING);


  GLfloat specularMat[] = {0.0, 0.0, 0.1, 1.0};
  GLfloat emissionMat[] = {0.0, 0.0, 0.0, 1.0};

  glMaterialfv(GL_FRONT, GL_SPECULAR, specularMat);
  glMaterialfv(GL_FRONT, GL_EMISSION, emissionMat);

  /*
  {
    glEnable(GL_LIGHT0);

    GLfloat ambientLight[] = {0., 0., 0., 1.0};
    GLfloat diffuseReflection[] = {0.8, 0.8, 0.8, 1.0};
    GLfloat specularLight[] = {0.5, 0.5, 0.5, 1.0};
    GLfloat pos[] = {  0.0, 0.0, -100.0, 1.0};
    GLfloat dir[] = {  0.0, 0.0,   1.0};

    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseReflection);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT0, GL_POSITION, pos);
    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, dir);
  }
    
  {
    glEnable(GL_LIGHT1);

    GLfloat ambientLight[] = {0.0, 0.0, 0.0, 1.0};
    GLfloat diffuseReflection[] = {0.2, 0.8, 0.2, 1.0};
    GLfloat specularLight[] = {0.2, 0.2, 0.2, 1.0};
    GLfloat pos[] = {-30., -30., 30.0, 1};
    GLfloat dir[] = {1/sqrt(3), 1/sqrt(3), -1/sqrt(3)};

    glLightfv(GL_LIGHT1, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuseReflection);
    glLightfv(GL_LIGHT1, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT1, GL_POSITION, pos);
    glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, dir);
  }

  */
 {
    glEnable(GL_LIGHT0);

    // (-30,   0,   0) - right
    // (  0, -30,   0) - up
    // (  0,   0, -30) - front
    Vec3 vpos(0.0, 0.0, -30.0);
    Vec3 vdir = vpos * (-1.0);
    vdir.normalize();

    GLfloat ambientLight[] = {0.0, 0.0, 0.0, 1.0};
    GLfloat diffuseReflection[] = {0.6, 0.6, 0.9, 1.0};
    GLfloat specularLight[] = {0.2, 0.2, 0.2, 1.0};
    GLfloat pos[] = {GLfloat(vpos.x()),  GLfloat(vpos.y()), GLfloat(vpos.z()), 1};
    GLfloat dir[] = {GLfloat(vdir.x()),  GLfloat(vdir.y()), GLfloat(vdir.z())};

    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseReflection);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT0, GL_POSITION, pos);
    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, dir);
  }

 

 {
    glEnable(GL_LIGHT1);

    // (-30,   0,   0) - right
    // (  0, -30,   0) - up
    // (  0,   0, -30) - front
    Vec3 vpos(-30.0, -20.0, 30.0);
    Vec3 vdir = vpos * (-1.0);
    vdir.normalize();

    GLfloat ambientLight[] = {0.0, 0.0, 0.0, 1.0};
    GLfloat diffuseReflection[] = {0.9, 0.6, 0.6, 1.0};
    GLfloat specularLight[] = {0.2, 0.2, 0.2, 1.0};
    GLfloat pos[] = {GLfloat(vpos.x()),  GLfloat(vpos.y()), GLfloat(vpos.z()), 1};
    GLfloat dir[] = {GLfloat(vdir.x()),  GLfloat(vdir.y()), GLfloat(vdir.z())};

    glLightfv(GL_LIGHT1, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuseReflection);
    glLightfv(GL_LIGHT1, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT1, GL_POSITION, pos);
    glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, dir);
  }

 
  {
    glEnable(GL_LIGHT2);

    // (-30,   0,   0) - right
    // (  0, -30,   0) - up
    // (  0,   0, -30) - front
    Vec3 vpos(30.0, -20.0, 30.0);
    Vec3 vdir = vpos * (-1.0);
    vdir.normalize();

    GLfloat ambientLight[] = {0.0, 0.0, 0.0, 1.0};
    GLfloat diffuseReflection[] = {0.6, 0.6, 0.8, 1.0};
    GLfloat specularLight[] = {0.2, 0.2, 0.2, 1.0};
    GLfloat pos[] = {GLfloat(vpos.x()),  GLfloat(vpos.y()), GLfloat(vpos.z()), 1};
    GLfloat dir[] = {GLfloat(vdir.x()),  GLfloat(vdir.y()), GLfloat(vdir.z())};

    glLightfv(GL_LIGHT2, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT2, GL_DIFFUSE, diffuseReflection);
    glLightfv(GL_LIGHT2, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT2, GL_POSITION, pos);
    glLightfv(GL_LIGHT2, GL_SPOT_DIRECTION, dir);
  }
 

}

void
GlView::resizeGL(int w, int h){
  int s = qMin(w, h);
  glViewport((w - h) / 2, (h - s) / 2, s, s);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glFrustum(-1.0, 1.0, -1.0, 1.0, 4, 12.0);
  glMatrixMode(GL_MODELVIEW);
}

double getDist(const Vector3d& n, int mode){
  MaterialTensor tensor = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  PolyMatrix pm = makeChristoffelPolyMatrix(tensor, n.x(), n.y(), n.z());
  Poly christ_poly = det(pm);
  Poly::RootVec r = christ_poly.all_roots();

  vector<double> ret(3);
  for(int t = 0; t < 3; ++t)
    ret[t] = sqrt(r[t].real() / rho);

  sort(ret.begin(), ret.end());
  
  for(int t = 0; t < 3; ++t)
    ret[t] = 1000. / ret[t]; 

  return ret[mode];
}

void
GlView::paintGL(){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glLoadIdentity();

  glTranslated(0, 0, -7);

  glRotated(xang, 1.0, 0.0, 0.0);
  glRotated(yang, 0.0, 1.0, 0.0);
  glRotated(zang, 0.0, 0.0, 1.0);

  double x;
  double y;
  double z;

  glPushAttrib(GL_CURRENT_BIT);
  glBegin(GL_LINES);{
    //Coordinate system
    glNormal3d(0,0,1);
    glColor3f(3.0, .0, 0.0);
    glVertex3f(-0.3, 0.0, 0.0);
    glVertex3f( 1.0, 0.0, 0.0);

    glNormal3d(1.0, 0, 0);
    glColor3f(0.0, 3.0, 0.0);
    glVertex3f(0.0, -0.3, 0.0);
    glVertex3f(0.0,  1.0, 0.0);

    glNormal3d(0, 1, 0);
    glColor3f(0.0, 0.0, 3.0);
    glVertex3f( 0.0, 0.0, -0.3);
    glVertex3f( 0.0, 0.0,  1.0);


    //Theta-phi radius from begin of coordinates
    double radTheta = theta / 360. * 2 * M_PI;
    double radPhi = phi / 360. * 2 * M_PI;
    x = sin(radTheta) * cos(radPhi);
    y = sin(radTheta) * sin(radPhi);
    z = cos(radTheta);
    glColor3f(3.0, 3.0, 0.0);
    glVertex3f( 0.0, 0.0, 0.0);
    glVertex3f( x, y,  z);

  }glEnd();
  glPopAttrib();

  

  //cylinder(Vector3d(x, y, z) * 0.6, Vector3d(x, y, z), 0.15, 0.01);
  //cone(Vector3d(x,y,z)*0.7, Vector3d(x,y,z), 0.05, 0.2);

  Vector3d vecWaveNorm(x,y,z);
  double dist = getDist(vecWaveNorm, pMainDialog -> getActiveSurface());
  yellowBall = vecWaveNorm * dist;
  draw_sphere(yellowBall.x(), yellowBall.y(), yellowBall.z(), 0.03, 5.0, 3.0, 0.0);

  cylinder(Vector3d(shiftX, shiftY, 0), Vector3d(0, 0, 1), 1, 0.005);

  
  double height = pHeightCalculator -> getHeight(shiftX / 1000, 
					     shiftY / 1000, 
					     pMainDialog -> getActiveSurface()) * 1000;
 
  pMainDialog->heightOut->setText(QString("%1").arg(height));
  draw_sphere(shiftX, 
	      shiftY, 
	      height, 
	      0.03, 
	      0, 
	      3, 
	      3);
  
  pWorld -> paint(isSurface);
  
}


void
GlView::changeXang(int xang_){
  xang = xang_;
  updateGL();
}

void
GlView::changeYang(int yang_){
  yang = yang_;
  updateGL();
}

void
GlView::changeZang(int zang_){
  zang = zang_;
  updateGL();
}

void
GlView::setVisible0(bool v){
  visible0 = v;
  updateGL();
}

void
GlView::setVisible1(bool v){
  visible1 = v;
  updateGL();
}

void
GlView::setVisible2(bool v){
  visible2 = v;
  updateGL();
}

void
GlView::changeTheta(int theta_){
  theta = theta_;
  pMainDialog->thetaVal->setText(QString("%1").arg(theta));
  updateGL();
}

void
GlView::changePhi(int phi_){
  phi = phi_;
  pMainDialog->phiVal->setText(QString("%1").arg(phi));
  updateGL();
}


void
GlView::setShiftx(int xshiftX){
  shiftX = double(xshiftX - 250) / 250;
  pMainDialog->shiftXVal->setText(QString("%1").arg(shiftX));
  updateGL();
}

void
GlView::setShifty(int xshiftY){
  shiftY = double(xshiftY - 250) / 250;
  pMainDialog->shiftYVal->setText(QString("%1").arg(shiftY));
  updateGL();
}

void
GlView::doTest(){
  cerr << "yellow ball position is " << yellowBall.x() << " " << yellowBall.y() << " " << yellowBall.z() << endl;
}

void GlView::showSlice() {
	QString stepStr = pMainDialog->sliceStepVal->text();
	vector<pair<double, double> > result = 
		pHeightCalculator->getAllHeights(shiftX / 1000., shiftY / 1000., 
		stepStr.toDouble() /1000.);
	//QMessageBox::information(pMainDialog, "Show Slice", 
	//	"Show slice with step "+stepStr);
	pGraph->setBody(result);
	pGraph->resize(320,200);
	pGraph->show();
	cout << "step=" << stepStr.toDouble() << endl;
}

void GlView::polarSlice() {
	QMessageBox::information(pMainDialog, "Polar Slice",
		"Polar slice with deg step "+pMainDialog->degStepVal->text());
	double degStep = pMainDialog->degStepVal->text().toDouble();
	if (degStep < 0.01) {
		 QMessageBox::information(pMainDialog, "Polar Slice",
		"deg step " + pMainDialog->degStepVal->text() + 
		" is too small or incorrect");
		return;
	}
	
    double radTheta = theta / 360. * 2 * M_PI;
    double radPhi = phi / 360. * 2 * M_PI;
    double x = sin(radTheta) * cos(radPhi);
    double y = sin(radTheta) * sin(radPhi);
    double z = cos(radTheta);

	Vec3 n(x,y,z);
	Vec3 v(1.,0.,0.);
	if ((theta > 45 && theta < 135) ||
		(theta > 225 && theta < 315)) {
		v = Vec3(0.,0.,1.);
	}
	Vec3 n1 = v & n;
	n1.normalize();
	Vec3 n2 = n & n1;
	n2.normalize(); //however, n2 already have norm 1 :)

	cout << "polarSlice:" << endl;
	cout << "n=" << n << endl;
	cout << "n1=" << n1  << endl;
	cout << "n2=" << n2 << endl;

	double radStep = degStep / 360. * 2 * M_PI;
  	MaterialTensor tensor = makeTetragonalMaterialTensor(
		5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
	vector<pair<double, double> > result;	

	for (double currAngle = 0.; currAngle < 2 * M_PI; currAngle += radStep) {
		Vec3 currN = n1*cos(currAngle)+n2*sin(currAngle);
		Poly::RootVec gammas = 
			calcGammas(det(makeChristoffelPolyMatrix(tensor, currN)));
		for (unsigned int i = 0; i < gammas.size(); ++i) {
			double gamma = real(gammas[i]); 
			double s = sqrt(rho/gamma);
			result.push_back(make_pair(cos(currAngle)*s, sin(currAngle)*s));
		}
	}
	pGraph->setBody(result);
	pGraph->resize(320,200);
	pGraph->show();
}

void
GlView::setSurfaceWire(int val){
  isSurface = 0 != val;
  updateGL();
}
