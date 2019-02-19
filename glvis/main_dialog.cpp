#include "main_dialog.h"

#include <QGLWidget>
#include "slow_gl_view.h"

#include <iostream>

using namespace std;

MainDialog::MainDialog(QWidget *parent, Qt::WindowFlags f) : QDialog(parent, f){
  setupUi(this);

  pglwidget = new SlowGlView;
  glLayout->addWidget(pglwidget);

  zdistSpinBox->setValue(pglwidget->zdist);

  connect(closeButton, SIGNAL(clicked()),
	  this, SLOT(close()));

  connect(zdistSpinBox, SIGNAL(valueChanged(double)),
	  pglwidget, SLOT(zdistChange(double)));

  connect(rotxSpinBox, SIGNAL(valueChanged(double)),
	  pglwidget, SLOT(xangleChange(double)));

  connect(rotySpinBox, SIGNAL(valueChanged(double)),
	  pglwidget, SLOT(yangleChange(double)));

  connect(rotzSpinBox, SIGNAL(valueChanged(double)),
	  pglwidget, SLOT(zangleChange(double)));

}
