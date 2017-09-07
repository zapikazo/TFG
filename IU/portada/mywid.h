#ifndef MYWID_H
#define MYWID_H

#include <QWidget>
#include <QPainter>
#include "experiment.h"

namespace Ui {
class MyWid;
}

class MyWid : public QWidget
{
    Q_OBJECT

public:
    explicit MyWid(QWidget *parent = 0);
    ~MyWid();

    std::vector<Experiment> experiments;
    std::vector<QColor> colors = {QColor("#FFFFFF"), QColor("#66B2FF"), QColor("#990099"), QColor("#00CC66"), QColor("#CC0000")};
    bool update;

    void updateWidget(){
        update = true;
    };
    void setExperiments(std::vector<Experiment> experiments);
    void paintEvent(QPaintEvent *event);

private:
    Ui::MyWid *ui;
};

#endif // MYWID_H
