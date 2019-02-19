#ifndef SLOW_GL_VIEW
#define SLOW_GL_VIEW

#include <QGLWidget>

class SlowGlView : public QGLWidget{
  Q_OBJECT
  public:
  double theta_;
  double phi_;
  double delta_;
  double zdist;


  SlowGlView(QWidget * = 0);

  void initializeGL();
  void resizeGL(int, int);
  void paintGL();


 public slots:
  void zdistChange(double);
  void xangleChange(double);
  void yangleChange(double);
  void zangleChange(double);
};

#endif
