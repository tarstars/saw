#include <iostream>

#include <QClipboard>
#include <QFileDialog>
#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>

#include "main_dialog.h"

using namespace std;

MainDialog::MainDialog(QWidget* parent, Qt::WindowFlags flags) : QDialog(parent, flags),
								 shx(50),
								 shy(50)
{
  setupUi(this);

  leaperture -> setText("0.002");

  connect(leaperture, SIGNAL(returnPressed()), this, SLOT(getText()));
  connect(pbLoadPicture, SIGNAL(clicked()), this, SLOT(openFile()));

  //setMouseTracking(true);
}

void
MainDialog::getText() {
  QString str = leaperture -> text();
  cout << str.toDouble() << endl;
}

void
MainDialog::openFile() {
  QString str = QFileDialog::getOpenFileName(this, tr("Open Image"), ".", tr("Image files (*.png *.jpg *.bmp)"));
  img.load(str);
  QSize sz = img.size();
  resize(sz.width() + shx, sz.height() + shy);
}

void
MainDialog::paintEvent(QPaintEvent*) {
  QPainter dc(this);

  dc.drawImage(shx, shy, img);
}

void
MainDialog::mouseMoveEvent(QMouseEvent *) {
}

void
MainDialog::mousePressEvent(QMouseEvent* pme) {
  int x = pme -> x();
  int y = pme -> y();
  double apert = leaperture -> text().toDouble();

  QSize sz = img.size();

  int mx = sz.width();
  int my = sz.height();

  if (mx == 0 || my == 0) {
    return;
  }
  
  double sx =  (double(x - shx) / mx - 0.5) * (2 * apert);
  double sy = -(double(y - shy) / my - 0.5) * (2 * apert);

  QString res = QString::number(sx) + ", " + QString::number(sy);
  lesxsy -> setText(res);
  QApplication::clipboard() -> setText(res);
}
