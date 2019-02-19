#ifndef GRAPHWIDGET_H
#define GRAPHWIDGET_H
#include <QWidget>
#include <QSize>
#include <vector>
#include <utility>

using namespace std;

class GraphWidget : public QWidget {
	Q_OBJECT
	public:
		GraphWidget(QWidget* parent = 0);
		void setInsets(const QSize&);
		void setBody(const vector<pair<double, double> >&);
		void setProportional(bool);
	private:
		vector<pair<double, double> > body;
		double maxX;
		double minX;
		double maxY;
		double minY;
		QSize insets;
		bool isProportional;
	protected:
		virtual void paintEvent(QPaintEvent*);
};

#endif
