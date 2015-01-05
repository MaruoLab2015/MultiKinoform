#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "bmplabel.h"
#include "dft.h"

#include <QDebug>
#include <QtMath>

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define alloc_error_check(p) { \
    if ((p) == NULL) { \
        fprintf(stderr, "Allocation Failure!\n"); \
        exit(1); \
    } \
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    kinoformArray(DynArray2D<double>(1024,1024)),// alloc kinoform array
    intensityArray(DynArray2D<double>(1024,1024)),
    convertCount(0),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

//    ui->bmpSizeComboBox->setCurrentIndex(2);
    getUIWidgetValue();
    init();
}

void MainWindow::init()
{
    QPixmap pixmap(kinoformAreaSize, kinoformAreaSize);
    confirmBmpLabel_ = new BmpLabel();
    confirmBmpLabel_->setFixedSize(pixmap.size());
    confirmBmpLabel_->setPixmap(pixmap);
    confirmBmpLabel_->show();

    initKinoform();
}

void MainWindow::initKinoform()
{
    for(int i=0;i<kinoformArray.width();i++)
        for(int j=0;j<kinoformArray.height();j++)
            kinoformArray.setValue(i,j,0);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_showIntensityButton_clicked()
{

    getUIWidgetValue();

    int n1, n2, n;

    int *alloc_1d_int(int n1);
    void free_1d_int(int *i);
    void free_1d_double(double *d);
    double *alloc_1d_double(int n1);
    double **alloc_2d_double(int n1, int n2);
    void free_2d_double(double **dd);
    n1 = kinoformAreaSize;                     //２次元高速DFTの標本数x
    n2 = 2*kinoformAreaSize;                   //２次元高速DFTの標本数y(実数と虚数で２倍)
    E = alloc_2d_double(n1, n2);
    E0 = alloc_2d_double(n1, n2);
    n = MAX(n1, n2 / 2);                                               //ビット反転のための作業領域
    ip = alloc_1d_int(2 + (int) sqrt(n + 0.5));                        //ビット反転のための作業領域
    n = MAX(n1 / 2, n2 / 4) + MAX(n1, n2);
    W = alloc_1d_double(n);                                            //三角関数テーブル
    ip[0] = 0;
    laserIntensity(beamWaist_,0.01);                                   //物体面の電場強度profileA[i][j]の計算

    qDebug() << E[0][0];

            int nanNum=0;
    for (int i=0;i<kinoformAreaSize;i++){
        for (int j=0;j<kinoformAreaSize;j++){


            if (!kinoformArray.getValue(i,j))
                nanNum++;
            E[i][2*j]=sqrt(intensityArray.getValue(i,j)*cos(kinoformArray.getValue(i,j)));	   //物体面の電場aa(j1,j2)の実数部(電場強度×cos(位相))
            E[i][2*j+1]=sqrt(intensityArray.getValue(i,j)*sin(kinoformArray.getValue(i,j)));     //物体面の電場aa(j1,j2)の虚数(電場強度×sin(位相))
        }
    }

    qDebug() << E[0][0];
    //２次元高速逆DTF
    cdft2d(n1, n2, 1, E, ip, W);                                       //２次元高速DTF 物体面の電場e(j1,j2)→ホログラム面の電場E(k1.k2)
    qDebug() << E[0][0];

    for (int i=0;i<kinoformAreaSize;i++){                                                      //行列の入れ替え
        for (int j=0;j<2*kinoformAreaSize;j++){
            E0[i][j]=E[i][j];
        }
    }

    for (int i=0;i<kinoformAreaSize;i++){                                                      //行列の入れ替え
        for (int j=0;j<2*kinoformAreaSize;j++){
            if(i<kinoformAreaSize/2 && j<kinoformAreaSize){
                E[i][j]=E0[i+kinoformAreaSize/2][j+kinoformAreaSize];
            }
            else if(i>=kinoformAreaSize/2 && j<kinoformAreaSize){
                E[i][j]=E0[i-kinoformAreaSize/2][j+kinoformAreaSize];
            }
            else if(i<kinoformAreaSize/2 && j>=kinoformAreaSize){
                E[i][j]=E0[i+kinoformAreaSize/2][j-kinoformAreaSize];
            }
            else{
                E[i][j]=E0[i-kinoformAreaSize/2][j-kinoformAreaSize];
            }
        }
    }

    for (int i=0;i<kinoformAreaSize;i++){
        for (int j=0;j<kinoformAreaSize;j++){
            intensityArray.setValue(i,j, E[i][2*j]*E[i][2*j]+E[i][2*j+1]*E[i][2*j+1]);//ホログラム面の電場強度b(k1,k2)=電場AA(k1,k2)の実数部^2+電場AA(k1,k2)の虚数部^2
        }
    }
    outputIntensity();                                                //ホログラム面の電場強度の表示
}

void MainWindow::on_showKinoformButton_clicked()
{
    getUIWidgetValue();
}

void MainWindow::on_convertButton_clicked()
{
    getUIWidgetValue();
    updateUIWidget();

    ui->numNumeratorLineEdit->setText(QString::number(convertCount+1));
    if (convertCount == 0)
        initKinoform();

    appliedLGBeam();
    appliedFlesnel();
    appliedFocalShift();
    outputKinoform();

    convertCount++;
    if ( convertCount == focusNumber)
    {
        convertCount = 0;
    }
}

void MainWindow::on_saveButton_clicked()
{

}


void MainWindow::getUIWidgetValue()
{
    waveLength_      = ui->waveLengthLineEdit->text().toDouble();
    focalLength_     = ui->focalLengthLineEdit->text().toDouble();
    beamWaist_       = ui->beamWaistLineEdit->text().toDouble();
    focusNumber      = ui->numLineEdit->text().toInt();
    tc_              = ui->tcLineEdit->text().toInt();
    initialPhase_l_  = ui->initialPhaseLLineEdit->text().toDouble();
    initialPhase_r_  = ui->initialPhaseRLineEdit->text().toDouble();
    initialPhase_z_  = ui->initialPhaseZLineEdit->text().toDouble();
    r_               = ui->rLineEdit->text().toDouble();
    phi_             = ui->phiLineEdit->text().toDouble();
    z_               = ui->zLineEdit->text().toDouble();
    kinoformAreaSize = ui->bmpSizeComboBox->currentText().toInt();
}

void MainWindow::updateUIWidget()
{
    confirmBmpLabel_->setFixedSize(kinoformAreaSize, kinoformAreaSize);
    ui->numDenominatorLineEdit->setText(QString::number(focusNumber));
}

void MainWindow::outputIntensity()
{
    QImage intensityImage = QImage(kinoformAreaSize,
                             kinoformAreaSize,
                             QImage::Format_ARGB32);

    // normarization of kinoform
    double MaxValue = 0;
    for (int i=0;i<kinoformAreaSize;i++)
        for (int j=0;j<kinoformAreaSize;j++)
            MaxValue = (intensityArray.getValue(i,j)>MaxValue) ? intensityArray.getValue(i,j):MaxValue; //　find the MaxValue

    for (int i=0;i<kinoformAreaSize;i++)
        for (int j=0;j<kinoformAreaSize;j++){
            int gray = qRound(intensityArray.getValue(i,j) / MaxValue * 255.0);
            intensityImage.setPixel(i,j, qRgb(gray,gray,gray));
        }

    confirmBmpLabel_->setPixmap(QPixmap::fromImage(intensityImage));

}

void MainWindow::outputKinoform()
{
    QImage kinoformImage = QImage(kinoformAreaSize,
                             kinoformAreaSize,
                             QImage::Format_ARGB32);


    // normarization of kinoform
    double MaxValue = 0;
    for (int i=0;i<kinoformAreaSize;i++)
        for (int j=0;j<kinoformAreaSize;j++)
            MaxValue = (kinoformArray.getValue(i,j)>MaxValue) ? kinoformArray.getValue(i,j):MaxValue; //　find the MaxValue

    for (int i=0;i<kinoformAreaSize;i++)
        for (int j=0;j<kinoformAreaSize;j++){
            int gray = qRound(kinoformArray.getValue(i,j) / MaxValue * 255.0);
            kinoformImage.setPixel(i,j, qRgb(gray,gray,gray));
        }

    confirmBmpLabel_->setPixmap(QPixmap::fromImage(kinoformImage));

    // 画面サイズとを合わせる
}

/////核関数///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MainWindow::appliedLGBeam()
{
    double phase = 0;
    for (int i=0;i<kinoformAreaSize;i++){
        for (int j=0;j<kinoformAreaSize;j++){
            if((i+j)%focusNumber == convertCount)
                phase=tc_*atan2((double)(i-(kinoformAreaSize/2.0)),(double)(j-(kinoformAreaSize/2.0)))+initialPhase_l_*M_PI/180.0;
            else
                phase=0.0;

            if ((kinoformArray.getValue(i,j) + phase) > 2*M_PI)
            {
                double tmpPhase = fmod((kinoformArray.getValue(i,j) + phase),2*M_PI);
                kinoformArray.setValue(i, j, 2*M_PI+tmpPhase);
            }
            else if( (kinoformArray.getValue(i,j) + phase) < 0){
                double tmpPhase = fmod(-1*(kinoformArray.getValue(i,j) + phase),2*M_PI);
                kinoformArray.setValue(i, j, 2*M_PI-tmpPhase);
            }
            else
                kinoformArray.setValue(i,j, kinoformArray.getValue(i,j) + phase);
        }
    }
}

void MainWindow::appliedFlesnel()
{
    double phase;

    for(int i=0;i<kinoformAreaSize;i++){
        for(int j=0;j<kinoformAreaSize;j++){
            if((i+j)%focusNumber == convertCount){
                double radius = ((double)(i-(kinoformAreaSize/2))*(double)(i-(kinoformAreaSize/2))
                        +(double)(j-(kinoformAreaSize/2))*(double)(j-(kinoformAreaSize/2)));
                phase=-(M_PI*radius*z_)/(focalLength_*waveLength_*(focalLength_+z_))+initialPhase_z_*M_PI/180.0;
            }else{
                phase=0.0;
            }
            kinoformArray.setValue(i,j,kinoformArray.getValue(i,j)+phase);                               //既存のキノフォームにフレネルレンズの位相を重ね合わせ
            if(kinoformArray.getValue(i,j) > 2*M_PI)
                kinoformArray.setValue(i,j,fmod(kinoformArray.getValue(i,j),2*M_PI));
            else if(kinoformArray.getValue(i,j) < 0){
                kinoformArray.setValue(i,j, fmod(-1*kinoformArray.getValue(i,j),2*M_PI));
                kinoformArray.setValue(i,j, -1*kinoformArray.getValue(i,j)+2*M_PI);
            }
        }
    }
}

void MainWindow::appliedFocalShift()
{
    double phase;
    for (int i=0;i<kinoformAreaSize;i++){
        for (int j=0;j<kinoformAreaSize;j++){
            if((i+j)%focusNumber == convertCount)
                phase = -2*M_PI*r_/(waveLength_*focalLength_)*(cos(phi_*M_PI/180.0)*(double)(i-(kinoformAreaSize/2))-sin(phi_*M_PI/180.0)*(double)(j-(kinoformAreaSize/2)))+initialPhase_r_*M_PI/180.0;
            else
                phase=0.0;
            kinoformArray.setValue(i,j,kinoformArray.getValue(i,j)+phase);
            if(kinoformArray.getValue(i,j) > 2*M_PI)                                      //焦点移動のキノフォーム
                kinoformArray.setValue(i,j,fmod(kinoformArray.getValue(i,j),2*M_PI));
            else if(kinoformArray.getValue(i,j) < 0){
                kinoformArray.setValue(i,j,fmod(-1*kinoformArray.getValue(i,j),2*M_PI));
                kinoformArray.setValue(i,j,2*M_PI-1*kinoformArray.getValue(i,j));
            }
        }
    }
}

/////光強度関連///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//レーザービームの強度分布
void MainWindow::laserIntensity(double Beamwaist, double laser_power)
{
    double radius;
    for (int i=0;i<kinoformAreaSize;i++){
        for (int j=0;j<kinoformAreaSize;j++){
            radius=((double)(i-(kinoformAreaSize/2.0))*(double)(i-(kinoformAreaSize/2.0))+(double)(j-(kinoformAreaSize/2.0))*(double)(j-(kinoformAreaSize/2.0)))*(20.0*20.0)/(kinoformAreaSize*kinoformAreaSize);
            intensityArray.setValue(i,j,(2.0*laser_power)*exp((-2*radius)/(Beamwaist*Beamwaist))/(M_PI*Beamwaist*Beamwaist));
        }
    }
}


void MainWindow::cdft2d(int n1, int n2, int isgn, double **a, int *ip, double *w)
{
    qDebug() << "dft";
    void makewt(int nw, int *ip, double *w);                                   //初期化？
    void bitrv2col(int n1, int n, int *ip, double **a);
    void bitrv2row(int n, int n2, int *ip, double **a);
    void cftbcol(int n1, int n, double **a, double *w);
    void cftbrow(int n, int n2, double **a, double *w);
    void cftfcol(int n1, int n, double **a, double *w);
    void cftfrow(int n, int n2, double **a, double *w);
    int n;

    n = n1 << 1;
    if (n < n2) {
        n = n2;
    }
    if (n > (ip[0] << 2)) {
        makewt(n >> 2, ip, w);
    }
    if (n2 > 4) {
        bitrv2col(n1, n2, ip + 2, a);
    }
    if (n1 > 2) {
        bitrv2row(n1, n2, ip + 2, a);
    }
    if (isgn < 0) {
        cftfcol(n1, n2, a, w);
        cftfrow(n1, n2, a, w);
    } else {
        cftbcol(n1, n2, a, w);
        cftbrow(n1, n2, a, w);
    }
}

