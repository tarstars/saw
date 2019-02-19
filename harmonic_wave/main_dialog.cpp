#include <QTimer>

#include "main_dialog.h"

#include "surface_view.h"

MainDialog::MainDialog(QWidget* parent, Qt::WindowFlags f):
  QDialog(parent, f){

  setupUi(this);

  glLayout -> addWidget(pGlView = new SurfaceView(this));
  pTimer = new QTimer(this);
  pTimer -> start(25);
  connect(pTimer, SIGNAL(timeout()),
	  pGlView, SLOT(onTimer()));

  connect(dialTheta, SIGNAL(valueChanged(int)),
	  pGlView, SLOT(setTheta(int)));
  dialTheta -> setValue(10);

  connect(dialPhi, SIGNAL(valueChanged(int)),
	  pGlView, SLOT(setPhi(int)));
  dialPhi -> setValue(10);

  connect(dialRotx, SIGNAL(valueChanged(int)),
	  pGlView, SLOT(setRotx(int)));
  connect(dialRoty, SIGNAL(valueChanged(int)),
	  pGlView, SLOT(setRoty(int)));
  connect(dialRotz, SIGNAL(valueChanged(int)),
	  pGlView, SLOT(setRotz(int)));

  dialRotx -> setValue(310);
  dialRoty -> setValue(0);
  dialRotz -> setValue(250);

  connect(showAngles, SIGNAL(clicked()),
	  pGlView, SLOT(showAngles()));
  
}

MainDialog::~MainDialog(){
  if (pGlView) 
    delete pGlView;

  if (pTimer)
    delete pTimer;
}
