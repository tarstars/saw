#ifndef MAIN_DIALOG
#define MAIN_DIALOG

#include <memory>
#include <Eigen/Geometry>

#include <QDialog>

#include "ui_main_dialog.h"

class GlView;
class World;

class MainDialog : public QDialog, public Ui::MainDialog{
  Q_OBJECT
    public:
  MainDialog(QWidget * = 0, Qt::WindowFlags = 0);

  GlView *pGlView;
  std::auto_ptr<World> pWorld;
  int getActiveSurface();
  Eigen::Vector3d getDirection();
  
};

#endif
