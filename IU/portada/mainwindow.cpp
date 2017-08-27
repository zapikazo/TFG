#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QProcess>

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
    //if (this->fileName.at(0) != null){
    //Comprobar si fileName tiene valor
        QProcess *process = new QProcess();
        QString program = "./program";
        QString file = this->fileName;
        process->start(program, QStringList() << file);
        //process->execute("./program", this->fileName);
        this->window.show();
   // }
}
