#ifndef SURFACE_VIEW
#define SURFACE_VIEW

#include <QGLWidget>

#include <eigen3/Eigen/Geometry>

class SurfaceView : public QGLWidget{
  Q_OBJECT
    
    double theta;
  double phi;
  Eigen::Vector3d dir;

  int rotx;
  int roty;
  int rotz;

    public:

  SurfaceView(QWidget* = 0);

 protected:

  void initializeGL();
  void resizeGL(int, int);
  void paintGL();

  public slots:
  
  void onTimer();
  void setTheta(int);
  void setPhi(int);

  void setRotx(int);
  void setRoty(int);
  void setRotz(int);

  void showAngles();
};

#endif
