#ifndef SECONDWINDOW_H
#define SECONDWINDOW_H

#include <QDialog>
#include <QString>
#include <QPainter>
#include <string>
#include <iostream>
#include <sstream>
#include <bitset>
#include "experiment.h"

namespace Ui {
class SecondWindow;
}

class SecondWindow : public QDialog
{
    Q_OBJECT

public:
    explicit SecondWindow(QWidget *parent = 0);
    ~SecondWindow();

    std::vector<Experiment> experiments;
    std::vector<QColor> colors = {QColor("#FFFFFF"), QColor("#66B2FF"), QColor("#990099"), QColor("#00CC66"), QColor("#CC0000")};

    std::pair<int, int> getCoordFromAddress(std::string number);
    void setExperiments(std::vector<Experiment> experiments);
    void setTitle(QString fileName){
        this->setWindowTitle(fileName);
    };

private slots:
    void on_comboBox_activated(int index);

private:
    Ui::SecondWindow *ui;
};

#endif // SECONDWINDOW_H
