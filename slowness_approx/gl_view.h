#ifndef GL_VIEW
#define GL_VIEW

#include <Eigen/Geometry>
#include <QGLWidget>
#include <memory>

#include "GraphWidget.h"

class World;
class HeightCalculator;
class MainDialog;

class GlView : public QGLWidget{
  Q_OBJECT

    public:

  std::auto_ptr<HeightCalculator> pHeightCalculator;

  double xang;
  double yang;
  double zang;

  Eigen::Vector3d yellowBall;

  int theta;
  int phi;

  double shiftX;
  double shiftY;

  bool visible0;
  bool visible1;
  bool visible2;

  bool isSurface;

  World *pWorld;
  MainDialog *pMainDialog;
  GraphWidget* pGraph;

  GlView(MainDialog * = 0);

  void initializeGL();
  void resizeGL(int, int);
  void paintGL();

  public slots:
  void changeXang(int);
  void changeYang(int);
  void changeZang(int);

  void setShiftx(int);
  void setShifty(int);
  
  void setVisible0(bool);
  void setVisible1(bool);
  void setVisible2(bool);

  void changeTheta(int);
  void changePhi(int);

  void doTest();
  void showSlice();
  void polarSlice();
  void setSurfaceWire(int);

};

#endif
