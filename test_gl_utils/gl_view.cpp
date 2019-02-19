#include <iostream>

#include "gl_view.h"

using namespace std;

GlView::GlView(QWidget* parent):
  QGLWidget(parent){

}

void
GlView::initializeGL(){
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  GLfloat lightPosition [] = {0.5, 5.0, 7.0, 1.0};
  glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
}

void
GlView::resizeGL(int w, int h){
  int s = qMin(w, h);
  glViewport((w - s) / 2, (h - s) / 2, s, s);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-1, 1, -1, 1, 4.0, 15.0);

  glMatrixMode(GL_MODELVIEW);

}

void
GlView::paintGL(){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  glTranslated(0, 0, -10);

  glBegin(GL_TRIANGLES);{
    glVertex3d(0, 1, 0);
    glVertex3d(-0.5, -0.5, 0);
    glVertex3d(0.5, -0.5, 0);
  }glEnd();
}

