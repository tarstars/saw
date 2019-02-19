#ifndef GL_VIEW
#define GL_VIEW

#include <QGLWidget>

class GlView : public QGLWidget {
  Q_OBJECT

    public:

  GlView(QWidget* = 0);

  void initializeGL();
  void resizeGL(int,int);
  void paintGL();
};

#endif
