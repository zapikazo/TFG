/********************************************************************************
** Form generated from reading UI file 'mywid.ui'
**
** Created by: Qt User Interface Compiler version 5.7.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MYWID_H
#define UI_MYWID_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MyWid
{
public:

    void setupUi(QWidget *MyWid)
    {
        if (MyWid->objectName().isEmpty())
            MyWid->setObjectName(QStringLiteral("MyWid"));
        MyWid->resize(400, 300);

        retranslateUi(MyWid);

        QMetaObject::connectSlotsByName(MyWid);
    } // setupUi

    void retranslateUi(QWidget *MyWid)
    {
        MyWid->setWindowTitle(QApplication::translate("MyWid", "Form", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class MyWid: public Ui_MyWid {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MYWID_H
