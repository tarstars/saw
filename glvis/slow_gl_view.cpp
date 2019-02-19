#include "slow_gl_view.h"

#include <iostream>
#include <QtOpenGL>


using namespace std;

SlowGlView::SlowGlView(QWidget *parent) : QGLWidget(parent), theta_(0), phi_(0), delta_(0), zdist(-5){
 
}

void
SlowGlView::initializeGL(){
  glEnable(GL_DEPTH_TEST);
  //glEnable(GL_CULL_FACE);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_BLEND);

  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  GLfloat ambientLight[] = {0.3f, 0.3f, 0.3f, 1.0f};
  GLfloat diffuseReflection[] = {0.2f, 0.5f, 0.3f, 1.0f};
  GLfloat pos[] = {3.0, 3.0, 0.0};
  GLfloat dir[] = {0.0, 0.0, 1.0};
  
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseReflection);
  glLightfv(GL_LIGHT0, GL_POSITION, pos);
  glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, dir);
}

void
SlowGlView::resizeGL(int w, int h){
  int s = qMin(w, h);
  glViewport((w - s) / 2, (h - s) / 2, s, s);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  //glOrtho(-4.0, +4.0, +4.0, -4.0, 4.0, 15.0);
  glFrustum(-1.0, +1.0, +1.0, -1.0, 4.0, 15.0);
  glMatrixMode(GL_MODELVIEW);
}

void
SlowGlView::paintGL(){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


  glLoadIdentity();

  glTranslated(GLfloat(0.0), GLfloat(0.0), GLfloat(zdist));

  glRotated(theta_,  1.0, 0.0, 0.0 );
  glRotated(phi_,    0.0, 1.0, 0.0 );
  glRotated(delta_,  0.0, 0.0, 1.0 );

  glBegin(GL_TRIANGLES);{
    glNormal3f( GLfloat( 0.0), GLfloat(0.0), GLfloat( -1));

    glColor3f( 1.0f, 0.0f, 0.0f);
    glVertex3f( GLfloat(-0.5), GLfloat(0.5), GLfloat(0.0));

    glColor3f( 0.0f, 1.0f, 0.0f);
    glVertex3f( GLfloat( 0.5), GLfloat(0.5), GLfloat(0.0));

    glColor3f( 0.0f, 0.0f, 1.0f);
    glVertex3f( GLfloat( 0.0), GLfloat(-0.5), GLfloat(0.0));
  }glEnd();
}

void
SlowGlView::zdistChange(double val){
  zdist = val;
  updateGL();
}

void
SlowGlView::xangleChange(double phi){
  phi_ = phi;
  updateGL();
}

void
SlowGlView::yangleChange(double theta){
  theta_ = theta;
  updateGL();
}


void
SlowGlView::zangleChange(double delta){
  delta_ = delta;
  updateGL();
}
