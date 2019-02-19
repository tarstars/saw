#include "main_dialog.h"

#include "gl_view.h"
#include "world.h"

#include <iostream>

using namespace Eigen;
using namespace std;

MainDialog::MainDialog(QWidget *parent, Qt::WindowFlags f):
  QDialog(parent, f), pWorld(new World){

  setupUi(this);

  pGlView = new GlView(this);
  glLayout -> addWidget(pGlView);

  pWorld -> registerWorld(pGlView);
  pWorld -> pMainDialog = this;

  connect(dialRotx, SIGNAL(valueChanged(int)), pGlView, SLOT(changeXang(int)));
  connect(dialRoty, SIGNAL(valueChanged(int)), pGlView, SLOT(changeYang(int)));
  connect(dialRotz, SIGNAL(valueChanged(int)), pGlView, SLOT(changeZang(int)));

  connect(dialTheta, SIGNAL(valueChanged(int)), pGlView, SLOT(changeTheta(int)));
  connect(dialPhi, SIGNAL(valueChanged(int)), pGlView, SLOT(changePhi(int)));

  connect(dialRotx, SIGNAL(valueChanged(int)), sbx, SLOT(setValue(int)));
  connect(dialRoty, SIGNAL(valueChanged(int)), sby, SLOT(setValue(int)));
  connect(dialRotz, SIGNAL(valueChanged(int)), sbz, SLOT(setValue(int)));

  connect(rbSelectSurface1, SIGNAL(toggled(bool)), pGlView, SLOT(setVisible0(bool)));
  connect(rbSelectSurface2, SIGNAL(toggled(bool)), pGlView, SLOT(setVisible1(bool)));
  connect(rbSelectSurface3, SIGNAL(toggled(bool)), pGlView, SLOT(setVisible2(bool)));

  connect(dialShiftX, SIGNAL(valueChanged(int)), pGlView, SLOT(setShiftx(int)));
  connect(dialShiftY, SIGNAL(valueChanged(int)), pGlView, SLOT(setShifty(int)));

  connect(bbDoTest, SIGNAL(clicked()), pGlView, SLOT(doTest()));
  connect(bbShowSlice, SIGNAL(clicked()), pGlView, SLOT(showSlice()));
  connect(bbPolarSlice, SIGNAL(clicked()), pGlView, SLOT(polarSlice()));	

  connect(cbWireSurface, SIGNAL(stateChanged(int)), pGlView, SLOT(setSurfaceWire(int)));
  cbWireSurface -> setChecked(true);


  dialRotx -> setValue(300);
  dialRotz -> setValue(250);

  dialTheta -> setValue(90);
  dialPhi -> setValue(45);

  rbSelectSurface1 -> setChecked(true);

  dialShiftX -> setValue(250);
  dialShiftY -> setValue(250);
  heightOut->setText(QString("%1").arg(0.731));
}

int
MainDialog::getActiveSurface(){
  if (rbSelectSurface1 -> isChecked())
    return 0;
  if (rbSelectSurface2 -> isChecked())
    return 1;
  return 2;
}

Vector3d
MainDialog::getDirection(){
  double radTheta = pGlView -> theta / 360. * 2 * M_PI;
  double radPhi = pGlView -> phi / 360. * 2 * M_PI;
  double x = sin(radTheta) * cos(radPhi);
  double y = sin(radTheta) * sin(radPhi);
  double z = cos(radTheta);

  return Vector3d(x,y,z);
}


