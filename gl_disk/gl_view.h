#ifndef GL_VIEW
#define GL_VIEW

#include <QGLWidget>

class GlView : public QGLWidget {
  Q_OBJECT
    public:

  int theta;
  int phi;

  int rotx, roty, rotz;

  GlView(QWidget* = 0);

  void initializeGL();
  void resizeGL(int, int);
  void paintGL();

  public slots:
  void setPhi(int phi);
  void setTheta(int theta);

  void setRotx(int);
  void setRoty(int);
  void setRotz(int);
};

#endif
