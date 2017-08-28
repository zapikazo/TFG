#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QProcess>
#include <QGraphicsScene>
#include <QGraphicsItem>
#include <QGraphicsRectItem>
#include <QSplitter>
#include "view.h"

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
        rellenoVentana();
        this->window.show();
    }
}

void MainWindow::rellenoVentana()
{
    this->window.setWindowTitle("MAIN");
    this->window.resize(680, 450);
    QHBoxLayout* mainLayout = new QHBoxLayout();
    QVBoxLayout* listLayout = new QVBoxLayout();
   // listLayout->setSizeConstraint(300);
    QVBoxLayout* graphLayout = new QVBoxLayout();
    //graphLayout->SetMaximumSize;
    createMemory();

    QPushButton* button = new QPushButton();
    QPushButton* button2 = new QPushButton();

    listLayout->addWidget(button);
    graphLayout->addWidget(button2);

    mainLayout->addLayout(graphLayout);
    mainLayout->addSpacing(50);
    mainLayout->addLayout(listLayout);

    this->window.setLayout(mainLayout);

    /*QGraphicsScene *scene;
    QSplitter *h1Splitter = new QSplitter;

    QSplitter *vSplitter = new QSplitter;
    vSplitter->setOrientation(Qt::Vertical);
    vSplitter->addWidget(h1Splitter);

    QGraphicsView *view = new QGraphicsView();
    view->setScene(scene);
    h1Splitter->addWidget(view);

    graphLayout->addWidget(vSplitter);*/
}

void MainWindow::createMemory()
{

}
