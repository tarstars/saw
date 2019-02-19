#include <QtGui>
#include <algorithm>
#include "GraphWidget.h"

bool comp1(const pair<double, double>& p1,
            const pair<double, double>& p2) {
    return p1.first < p2.first;
}

bool comp2(const pair<double, double>& p1,
            const pair<double, double>& p2) {
    return p1.second < p2.second;
}

GraphWidget::GraphWidget(QWidget* parent): QWidget(parent),
	maxX(1.), minX(0.), maxY(1.), minY(0.), insets(20,20),
	isProportional(true) {
	QPalette palette;
	palette.setColor(backgroundRole(),
		Qt::black);
	setPalette(palette);
}

void GraphWidget::setInsets(const QSize& insets) {
	this->insets = insets;
}

void GraphWidget::setBody(const vector<pair<double, double> >& body) {
	this->body = body;
	maxX = max_element(body.begin(), body.end(), comp1)->first;	
	maxY = max_element(body.begin(), body.end(), comp2)->second;
	minX = min_element(body.begin(), body.end(), comp1)->first;
	minY = min_element(body.begin(), body.end(), comp2)->second;
	if (maxX == minX) maxX = maxX+1.;
	if (maxY == minY) maxY = maxY+1.;
}

void GraphWidget::setProportional(bool prop) {
	isProportional = prop;
}

void GraphWidget::paintEvent(QPaintEvent*) {
	QPainter painter(this);
	painter.setPen(Qt::white);

	int width = size().width()-2*insets.width();
	int height = size().height()-2*insets.height();
	double kx = width/(maxX-minX);
	double ky = height/(maxY-minY);
	if (isProportional) {
		kx = min(kx,ky);
		ky = kx;
	}
	for (vector<pair<double, double> >::iterator i = body.begin();
			i != body.end(); ++i) {
		
		double x = i->first;
		double y = i->second;

		int x1 = (int) ((x-minX)*kx + insets.width());
		int y1 = (int) (height - (y-minY)*ky + insets.height());

		painter.drawPoint(x1, y1);
	}
}
