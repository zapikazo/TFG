#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QProcess>
#include <fstream>
#include <iostream>
#include "view.h"
#include "secondwindow.h"
#include "mywid.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->setWindowTitle("SBUs vs. MCUs");
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    QStringList fileName = QFileDialog::getOpenFileNames(this, "Select a file to open...", QDir::currentPath());
    if(fileName.size() != 0){
        ui->label->clear();
        ui->label->setText(fileName.at(0));

        // Parsing
        QString route = fileName.at(0);
        QStringList pieces = route.split( "/" );
        // MainWindow.fileName has the correct name
        this->fileName = pieces.value( pieces.length() - 1 );
    }
}


void MainWindow::on_pushButton_2_clicked()
{
    if (this->fileName.size() != 0){  // We can start the program if fileName has some value
        QProcess *process = new QProcess();
        QString program = "./program";
        QString file = this->fileName;
        process->start(program, QStringList() << file);
        //process->execute("./program", this->fileName);

        SecondWindow secondWindow;
        secondWindow.setTitle(file);
        secondWindow.setExperiments(fileReading()); //Devuelve el vector de experimentos
        secondWindow.setModal(true);
        secondWindow.exec();
    }
}

std::vector<Experiment> MainWindow::fileReading()
{
    std::vector<Experiment> exp;
    int exps = 0;
    int maxMCU = 0;
    std::vector<std::pair<int, int>> events;

    std::ifstream ifs ("results.txt", std::ifstream::in);
    //Lectura de fichero
    ifs >> exps;
    ifs >> maxMCU;
    for (int i = 0; i < exps; i++){
        int n;
        ifs >> n;
        events.push_back(std::pair<int, int>(1, n));
        for (int j = 2; j <= maxMCU; j++){
            ifs >> n;
            ifs >> n;
            events.push_back(std::pair<int, int>(j, n));
        }
        Experiment e;
        e.events = events;
        exp.push_back(e);
        events.clear();
    }

    /*int32_t*** result;
    result = new int32_t**[this->maxRows];
    for(int i = 0; i < this->maxRows; i++) {
        result[i] = new int32_t*[4];
        for(int j = 0; j < 4; j++){
            result[i][j] = new int32_t[this->experiments];
        }
    }*/



    ifs.close();
    return exp;
}
