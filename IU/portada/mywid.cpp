#include "mywid.h"
#include "ui_mywid.h"

MyWid::MyWid(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::MyWid)
{
    ui->setupUi(this);
    update = true;
}

MyWid::~MyWid()
{
    delete ui;
}

/*QPen getPenColor(int color){
    switch(color){
        case 1:
            return QPen penG(Qt::green, 30, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
        case 2:
            return QPen penB(Qt::blue, 30, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
        case 3:
            return QPen penR(Qt::red, 30, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
        case 4:
            return QPen penBK(Qt::black, 30, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
        default:
            return QPen penBK(Qt::black, 30, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
    }
}
*/
void MyWid::paintEvent(QPaintEvent *event)
{
   /* if (!update) {
        return;
    }*/
    QPainter painter(this);

    // Variables de la memoria
    int nLineas = 100;		// Lineas de memoria
    int nBloques = 64;		// Bloques por linea
    int nIOs = 8;			// IOs por bloque
    int nDirecciones = 8;	// Direcciones por IO

    // Variables de la matriz
    int iX = 0;				// Indice x de la matriz
    int iY = 0;				// Indice y de la matriz

    // Variables de dibujo
    int origX = 50;			// Posicion x de mas a la izquierda
    int x = origX;			// Posicion x de la primera dir.
    int y = 50;				// Posicion y de la primera dir.
    int szX = 3;			// Ancho de los bits (direcciones)
    int szY = 3;			// Alto de los bits
    int bitGap = 2;			// Espacio entre bits	(Horiz)
    int ioGap = 4;			// Espacio entre IOS	(Horiz)
    int blockGap = 8;		// Espacio entre bloques (Horiz)
    int rowGap = 2;			// Espacio entre lineas (Vert)

    for(int i=0; i < nLineas; ++i){
        for(int j=0; j < nBloques; ++j){
            for(int k=0; k < nIOs; ++k){
                for(int l=0; l < nDirecciones; ++l){
                   // painter.setPen(getPenColor(matriz[iX][iY]);
                    painter.drawRect(x,y, szX,szY);
                    if(l != nDirecciones-1)
                        x += szX + bitGap; // El siguiente se pone 'bitGap' pixeles a la derecha de donde termina el anterior (szX)
                }
                if(k != nIOs-1)
                    x += szX + ioGap;	// El IO termina y se separa del siguiente
            }
            if(j != nBloques-1)
                x += szX + blockGap;	// El bloque termina y se separa del siguiente
        }
        y += szY + rowGap;	// Linea terminada, se baja a la switch
        x = origX;		// Se vuelve a empezar por la izquierda
    }

   /* for(int i = 0; i < 4096; i++){
        painter.drawRect(10, (i*10)+10, 4096, 1);
    }

    for(int i = 0; i < 4096; i++){
        painter.drawRect((i*10)+10, 10, 1, 4096);
    }*/

    painter.end();
    update = false;
}



