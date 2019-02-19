/********************************************************************************
** Form generated from reading UI file 'main_dialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.5
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAIN_DIALOG_H
#define UI_MAIN_DIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDial>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainDialog
{
public:
    QDialogButtonBox *dialogButtons;
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *glLayout;
    QDial *dialRotx;
    QLabel *label;
    QLabel *label_2;
    QDial *dialRoty;
    QLabel *label_3;
    QDial *dialRotz;
    QSpinBox *sbx;
    QSpinBox *sby;
    QSpinBox *sbz;
    QRadioButton *rbSelectSurface1;
    QRadioButton *rbSelectSurface2;
    QRadioButton *rbSelectSurface3;
    QDial *dialTheta;
    QDial *dialPhi;
    QLabel *label_4;
    QLabel *label_5;
    QPushButton *bbDoTest;
    QDial *dialShiftX;
    QDial *dialShiftY;
    QLabel *label_6;
    QLabel *label_7;
    QCheckBox *cbWireSurface;
    QLabel *label_8;
    QLineEdit *heightOut;
    QLineEdit *shiftXVal;
    QLineEdit *shiftYVal;
    QPushButton *bbShowSlice;
    QLabel *label_9;
    QLineEdit *sliceStepVal;
    QLineEdit *thetaVal;
    QLineEdit *phiVal;
    QLabel *label_10;
    QLineEdit *degStepVal;
    QPushButton *bbPolarSlice;
    QButtonGroup *bgSurfaceSelect;

    void setupUi(QDialog *MainDialog)
    {
        if (MainDialog->objectName().isEmpty())
            MainDialog->setObjectName(QStringLiteral("MainDialog"));
        MainDialog->resize(1007, 680);
        dialogButtons = new QDialogButtonBox(MainDialog);
        dialogButtons->setObjectName(QStringLiteral("dialogButtons"));
        dialogButtons->setGeometry(QRect(750, 480, 191, 31));
        dialogButtons->setOrientation(Qt::Horizontal);
        dialogButtons->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        horizontalLayoutWidget = new QWidget(MainDialog);
        horizontalLayoutWidget->setObjectName(QStringLiteral("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(30, 20, 631, 611));
        glLayout = new QHBoxLayout(horizontalLayoutWidget);
        glLayout->setObjectName(QStringLiteral("glLayout"));
        glLayout->setContentsMargins(0, 0, 0, 0);
        dialRotx = new QDial(MainDialog);
        dialRotx->setObjectName(QStringLiteral("dialRotx"));
        dialRotx->setGeometry(QRect(720, 30, 70, 70));
        dialRotx->setMaximum(360);
        dialRotx->setValue(0);
        label = new QLabel(MainDialog);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(740, 10, 41, 16));
        label_2 = new QLabel(MainDialog);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setGeometry(QRect(830, 10, 59, 15));
        dialRoty = new QDial(MainDialog);
        dialRoty->setObjectName(QStringLiteral("dialRoty"));
        dialRoty->setGeometry(QRect(810, 30, 70, 70));
        dialRoty->setMaximum(360);
        label_3 = new QLabel(MainDialog);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setGeometry(QRect(930, 10, 59, 15));
        dialRotz = new QDial(MainDialog);
        dialRotz->setObjectName(QStringLiteral("dialRotz"));
        dialRotz->setGeometry(QRect(910, 30, 70, 70));
        dialRotz->setMaximum(360);
        dialRotz->setValue(0);
        sbx = new QSpinBox(MainDialog);
        sbx->setObjectName(QStringLiteral("sbx"));
        sbx->setGeometry(QRect(730, 110, 54, 23));
        sbx->setMaximum(360);
        sby = new QSpinBox(MainDialog);
        sby->setObjectName(QStringLiteral("sby"));
        sby->setGeometry(QRect(820, 110, 54, 23));
        sby->setMaximum(360);
        sbz = new QSpinBox(MainDialog);
        sbz->setObjectName(QStringLiteral("sbz"));
        sbz->setGeometry(QRect(920, 110, 54, 23));
        sbz->setMaximum(360);
        rbSelectSurface1 = new QRadioButton(MainDialog);
        bgSurfaceSelect = new QButtonGroup(MainDialog);
        bgSurfaceSelect->setObjectName(QStringLiteral("bgSurfaceSelect"));
        bgSurfaceSelect->addButton(rbSelectSurface1);
        rbSelectSurface1->setObjectName(QStringLiteral("rbSelectSurface1"));
        rbSelectSurface1->setGeometry(QRect(730, 160, 93, 21));
        rbSelectSurface2 = new QRadioButton(MainDialog);
        bgSurfaceSelect->addButton(rbSelectSurface2);
        rbSelectSurface2->setObjectName(QStringLiteral("rbSelectSurface2"));
        rbSelectSurface2->setGeometry(QRect(730, 190, 93, 21));
        rbSelectSurface3 = new QRadioButton(MainDialog);
        bgSurfaceSelect->addButton(rbSelectSurface3);
        rbSelectSurface3->setObjectName(QStringLiteral("rbSelectSurface3"));
        rbSelectSurface3->setGeometry(QRect(730, 220, 93, 21));
        dialTheta = new QDial(MainDialog);
        dialTheta->setObjectName(QStringLiteral("dialTheta"));
        dialTheta->setGeometry(QRect(830, 180, 70, 70));
        dialTheta->setMaximum(360);
        dialTheta->setValue(0);
        dialTheta->setWrapping(true);
        dialPhi = new QDial(MainDialog);
        dialPhi->setObjectName(QStringLiteral("dialPhi"));
        dialPhi->setGeometry(QRect(920, 180, 70, 70));
        dialPhi->setMaximum(360);
        dialPhi->setWrapping(true);
        label_4 = new QLabel(MainDialog);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setGeometry(QRect(850, 160, 41, 16));
        label_5 = new QLabel(MainDialog);
        label_5->setObjectName(QStringLiteral("label_5"));
        label_5->setGeometry(QRect(950, 160, 31, 16));
        bbDoTest = new QPushButton(MainDialog);
        bbDoTest->setObjectName(QStringLiteral("bbDoTest"));
        bbDoTest->setGeometry(QRect(720, 290, 71, 22));
        dialShiftX = new QDial(MainDialog);
        dialShiftX->setObjectName(QStringLiteral("dialShiftX"));
        dialShiftX->setGeometry(QRect(830, 280, 70, 70));
        dialShiftX->setMaximum(500);
        dialShiftX->setValue(0);
        dialShiftY = new QDial(MainDialog);
        dialShiftY->setObjectName(QStringLiteral("dialShiftY"));
        dialShiftY->setGeometry(QRect(920, 280, 70, 70));
        dialShiftY->setMaximum(500);
        label_6 = new QLabel(MainDialog);
        label_6->setObjectName(QStringLiteral("label_6"));
        label_6->setGeometry(QRect(820, 280, 16, 16));
        label_7 = new QLabel(MainDialog);
        label_7->setObjectName(QStringLiteral("label_7"));
        label_7->setGeometry(QRect(920, 280, 16, 16));
        cbWireSurface = new QCheckBox(MainDialog);
        cbWireSurface->setObjectName(QStringLiteral("cbWireSurface"));
        cbWireSurface->setGeometry(QRect(730, 260, 70, 17));
        label_8 = new QLabel(MainDialog);
        label_8->setObjectName(QStringLiteral("label_8"));
        label_8->setGeometry(QRect(710, 440, 51, 17));
        heightOut = new QLineEdit(MainDialog);
        heightOut->setObjectName(QStringLiteral("heightOut"));
        heightOut->setEnabled(false);
        heightOut->setGeometry(QRect(760, 440, 201, 21));
        shiftXVal = new QLineEdit(MainDialog);
        shiftXVal->setObjectName(QStringLiteral("shiftXVal"));
        shiftXVal->setEnabled(false);
        shiftXVal->setGeometry(QRect(830, 350, 71, 27));
        shiftYVal = new QLineEdit(MainDialog);
        shiftYVal->setObjectName(QStringLiteral("shiftYVal"));
        shiftYVal->setEnabled(false);
        shiftYVal->setGeometry(QRect(920, 350, 71, 27));
        bbShowSlice = new QPushButton(MainDialog);
        bbShowSlice->setObjectName(QStringLiteral("bbShowSlice"));
        bbShowSlice->setGeometry(QRect(710, 380, 81, 21));
        label_9 = new QLabel(MainDialog);
        label_9->setObjectName(QStringLiteral("label_9"));
        label_9->setGeometry(QRect(800, 380, 41, 16));
        sliceStepVal = new QLineEdit(MainDialog);
        sliceStepVal->setObjectName(QStringLiteral("sliceStepVal"));
        sliceStepVal->setGeometry(QRect(840, 380, 121, 21));
        thetaVal = new QLineEdit(MainDialog);
        thetaVal->setObjectName(QStringLiteral("thetaVal"));
        thetaVal->setEnabled(false);
        thetaVal->setGeometry(QRect(840, 250, 51, 21));
        phiVal = new QLineEdit(MainDialog);
        phiVal->setObjectName(QStringLiteral("phiVal"));
        phiVal->setEnabled(false);
        phiVal->setGeometry(QRect(930, 250, 51, 21));
        label_10 = new QLabel(MainDialog);
        label_10->setObjectName(QStringLiteral("label_10"));
        label_10->setGeometry(QRect(800, 410, 61, 16));
        degStepVal = new QLineEdit(MainDialog);
        degStepVal->setObjectName(QStringLiteral("degStepVal"));
        degStepVal->setGeometry(QRect(870, 410, 91, 21));
        bbPolarSlice = new QPushButton(MainDialog);
        bbPolarSlice->setObjectName(QStringLiteral("bbPolarSlice"));
        bbPolarSlice->setGeometry(QRect(710, 410, 81, 21));

        retranslateUi(MainDialog);
        QObject::connect(dialogButtons, SIGNAL(accepted()), MainDialog, SLOT(accept()));
        QObject::connect(dialogButtons, SIGNAL(rejected()), MainDialog, SLOT(reject()));

        QMetaObject::connectSlotsByName(MainDialog);
    } // setupUi

    void retranslateUi(QDialog *MainDialog)
    {
        MainDialog->setWindowTitle(QApplication::translate("MainDialog", "Dialog", Q_NULLPTR));
        label->setText(QApplication::translate("MainDialog", "RotX", Q_NULLPTR));
        label_2->setText(QApplication::translate("MainDialog", "RotY", Q_NULLPTR));
        label_3->setText(QApplication::translate("MainDialog", "RotZ", Q_NULLPTR));
        rbSelectSurface1->setText(QApplication::translate("MainDialog", "Surface1", Q_NULLPTR));
        rbSelectSurface2->setText(QApplication::translate("MainDialog", "Surface2", Q_NULLPTR));
        rbSelectSurface3->setText(QApplication::translate("MainDialog", "Surface3", Q_NULLPTR));
        label_4->setText(QApplication::translate("MainDialog", "Theta", Q_NULLPTR));
        label_5->setText(QApplication::translate("MainDialog", "Phi", Q_NULLPTR));
        bbDoTest->setText(QApplication::translate("MainDialog", "Do test", Q_NULLPTR));
        label_6->setText(QApplication::translate("MainDialog", "X", Q_NULLPTR));
        label_7->setText(QApplication::translate("MainDialog", "Y", Q_NULLPTR));
        cbWireSurface->setText(QApplication::translate("MainDialog", "Surface", Q_NULLPTR));
        label_8->setText(QApplication::translate("MainDialog", "height:", Q_NULLPTR));
        bbShowSlice->setText(QApplication::translate("MainDialog", "Show Slice", Q_NULLPTR));
        label_9->setText(QApplication::translate("MainDialog", "step:", Q_NULLPTR));
        label_10->setText(QApplication::translate("MainDialog", "deg step:", Q_NULLPTR));
        bbPolarSlice->setText(QApplication::translate("MainDialog", "Polar Slice", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class MainDialog: public Ui_MainDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAIN_DIALOG_H
