#include "secondwindow.h"
#include "ui_secondwindow.h"

SecondWindow::SecondWindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SecondWindow)
{
    ui->setupUi(this);
}

SecondWindow::~SecondWindow()
{
    delete ui;
}

void SecondWindow::setExperiments(std::vector<Experiment> experiments)
{
    ui->comboBox->addItem("All tests");
    this->experiments = experiments;
    for (size_t i = 0; i < experiments.size(); i++) {
        QString s;
        s.setNum(i + 1);
        s.prepend("Test ");
        ui->comboBox->addItem(s);
    }

    ui->tableWidget->resizeColumnsToContents();
    ui->tableWidget->setColumnWidth(2, 140);
    ui->tableWidget->verticalHeader()->setVisible(false);

    ui->tableWidget->setRowCount(experiments[0].events.size());
    for (size_t i = 0; i <  experiments[0].events.size(); i++) {
        QTableWidgetItem* item = new QTableWidgetItem(" ");
        item->setBackgroundColor(colors[i]);
        ui->tableWidget->setItem(i, 0, item);
        QString s;
        s.setNum(i + 1);
        s.append("-bit");
        ui->tableWidget->setItem(i, 1, new QTableWidgetItem(s));
        int count = 0;
        for (size_t j = 0; j < experiments.size(); j++) {
            count += experiments[j].events[i].second;
        }
        s.setNum(count);
        ui->tableWidget->setItem(i, 2, new QTableWidgetItem(s));
    }

}

std::pair<int, int> SecondWindow::getCoordFromAddress(std::string number)
{
    std::stringstream ss;
    ss << std::hex << number;
    unsigned n;
    ss >> n;
    std::bitset<24> b(n);
    std::cout << b.to_string() << std::endl;
    std::cout << "sapin" << b[0] << std::endl;
    std::cout << "ax0 ax10" << b[0] << std::endl;

    std::cout << "sapin" << b[0] << std::endl;

    std::cout << "sapin" << b[0] << std::endl;


    //Colocacion de cypress


    //ultima posición zona memoria y Bin a dec primeras 11 posiciones

    //12 a 20 columna

}


void SecondWindow::on_comboBox_activated(int index)
{
    if(index == 0){
        for (size_t i = 0; i <  experiments[0].events.size(); i++) {
            QTableWidgetItem* item = new QTableWidgetItem(" ");
            item->setBackgroundColor(colors[i]);
            ui->tableWidget->setItem(i, 0, item);
            QString s;
            s.setNum(i + 1);
            s.append("-bit");
            ui->tableWidget->setItem(i, 1, new QTableWidgetItem(s));
            int count = 0;
            for (size_t j = 0; j < experiments.size(); j++) {
                count += experiments[j].events[i].second;
            }
            s.setNum(count);
            ui->tableWidget->setItem(i, 2, new QTableWidgetItem(s));
        }
    } else {
        for (size_t i = 0; i <  experiments[0].events.size(); i++) {
            QTableWidgetItem* item = new QTableWidgetItem(" ");
            item->setBackgroundColor(colors[i]);
            ui->tableWidget->setItem(i, 0, item);
            QString s;
            s.setNum(i + 1);
            s.append("-bit");
            ui->tableWidget->setItem(i, 1, new QTableWidgetItem(s));
            s.setNum(experiments[index-1].events[i].second);
            ui->tableWidget->setItem(i, 2, new QTableWidgetItem(s));
        }
    }
}

