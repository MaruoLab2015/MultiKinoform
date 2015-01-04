#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "bmplabel.h"

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

    ui->bmpSizeComboBox->setCurrentIndex(2);
    getUIWidgetValue();
    init();
}

void MainWindow::init()
{
//    QPixmap pixmap("/Users/genki/holo.bmp");
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

    for (int i=0;i<kinoformAreaSize;i++){
        for (int j=0;j<kinoformAreaSize;j++){
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

    ui->numNumeratorLineEdit->setText(QString::number(convertCount));
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
    QImage intensityImage = QImage(intensityArray.width(),
                             intensityArray.height(),
                             QImage::Format_ARGB32);

    // normarization of kinoform
    double MaxValue = 0;
    for (int i=0;i<intensityArray.width();i++)
        for (int j=0;j<intensityArray.height();j++)
            MaxValue = (intensityArray.getValue(i,j)>MaxValue) ? intensityArray.getValue(i,j):MaxValue; //　find the MaxValue

    for (int i=0;i<intensityArray.width();i++)
        for (int j=0;j<intensityArray.height();j++){
            int gray = qRound(intensityArray.getValue(i,j) / MaxValue * 255.0);
            intensityImage.setPixel(i,j, qRgb(gray,gray,gray));
        }

    confirmBmpLabel_->setPixmap(QPixmap::fromImage(intensityImage));

}

void MainWindow::outputKinoform()
{
    QImage kinoformImage = QImage(kinoformArray.width(),
                             kinoformArray.height(),
                             QImage::Format_ARGB32);

    // normarization of kinoform
    double MaxValue = 0;
    for (int i=0;i<kinoformArray.width();i++)
        for (int j=0;j<kinoformArray.height();j++)
            MaxValue = (kinoformArray.getValue(i,j)>MaxValue) ? kinoformArray.getValue(i,j):MaxValue; //　find the MaxValue

    for (int i=0;i<kinoformArray.width();i++)
        for (int j=0;j<kinoformArray.height();j++){
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
    for (int i=0;i<kinoformArray.width();i++){
        for (int j=0;j<kinoformArray.height();j++){
            if((i+j)%focusNumber == convertCount)
                phase=tc_*atan2((double)(i-(kinoformArray.width()/2)),(double)(j-(kinoformArray.height()/2)))+initialPhase_l_*M_PI/180.0;
            else
                phase=0.0;

            if ((kinoformArray.getValue(i,j) + phase) > 2*M_PI)
                kinoformArray.setValue(i, j, fmod(kinoformArray.getValue(i, j), 2*M_PI));
            else if( (kinoformArray.getValue(i,j) + phase) < 0){
                double tmpPhase = fmod(-1*(kinoformArray.getValue(i,j) + phase),2*M_PI);
                kinoformArray.setValue(i, j, 2*M_PI-1*tmpPhase);
            }
            else
                kinoformArray.setValue(i,j, kinoformArray.getValue(i,j) + phase);
        }
    }
}

void MainWindow::appliedFlesnel()
{
    double phase;

    for(int i=0;i<kinoformArray.width();i++){
        for(int j=0;j<kinoformArray.height();j++){
            if((i+j)%focusNumber == convertCount){
                double radius = ((double)(i-(kinoformArray.width()/2))*(double)(i-(kinoformArray.width()/2))
                        +(double)(j-(kinoformArray.height()/2))*(double)(j-(kinoformArray.height()/2)));
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
    for (int i=0;i<kinoformArray.width();i++){
        for (int j=0;j<kinoformArray.height();j++){
            if((i+j)%focusNumber == convertCount)
                phase = -2*M_PI*r_/(waveLength_*focalLength_)*(cos(phi_*M_PI/180.0)*(double)(i-(kinoformArray.width()/2))-sin(phi_*M_PI/180.0)*(double)(j-(kinoformArray.height()/2)))+initialPhase_r_*M_PI/180.0;
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
    for (int i=0;i<intensityArray.width();i++){
        for (int j=0;j<intensityArray.height();j++){
            radius=((double)(i-(intensityArray.width()/2))*(double)(i-(intensityArray.width()/2))+(double)(j-(intensityArray.height()/2))*(double)(j-(intensityArray.height()/2)))*(20.0*20.0)/(intensityArray.width()*intensityArray.height());
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

void makewt(int nw, int *ip, double *w)
{
    void bitrv2(int n, int *ip, double *a);
    int nwh, j;
    double delta, x, y;

    ip[0] = nw;
    ip[1] = 1;
    if (nw > 2) {
        nwh = nw >> 1;
        delta = atan(1.0) / nwh;
        w[0] = 1;
        w[1] = 0;
        w[nwh] = cos(delta * nwh);
        w[nwh + 1] = w[nwh];
        for (j = 2; j <= nwh - 2; j += 2) {
            x = cos(delta * j);
            y = sin(delta * j);
            w[j] = x;
            w[j + 1] = y;
            w[nw - j] = y;
            w[nw - j + 1] = x;
        }
        bitrv2(nw, ip + 2, w);
    }
}

/* -------- child routines -------- */

void bitrv2(int n, int *ip, double *a)
{
    int j, j1, k, k1, l, m, m2;
    double xr, xi;

    ip[0] = 0;
    l = n;
    m = 1;
    while ((m << 2) < l) {
        l >>= 1;
        for (j = 0; j <= m - 1; j++) {
            ip[m + j] = ip[j] + l;
        }
        m <<= 1;
    }
    if ((m << 2) > l) {
        for (k = 1; k <= m - 1; k++) {
            for (j = 0; j <= k - 1; j++) {
                j1 = (j << 1) + ip[k];
                k1 = (k << 1) + ip[j];
                xr = a[j1];
                xi = a[j1 + 1];
                a[j1] = a[k1];
                a[j1 + 1] = a[k1 + 1];
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
        }
    } else {
        m2 = m << 1;
        for (k = 1; k <= m - 1; k++) {
            for (j = 0; j <= k - 1; j++) {
                j1 = (j << 1) + ip[k];
                k1 = (k << 1) + ip[j];
                xr = a[j1];
                xi = a[j1 + 1];
                a[j1] = a[k1];
                a[j1 + 1] = a[k1 + 1];
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += m2;
                xr = a[j1];
                xi = a[j1 + 1];
                a[j1] = a[k1];
                a[j1 + 1] = a[k1 + 1];
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
        }
    }
}


void bitrv2col(int n1, int n, int *ip, double **a)
{
    int i, j, j1, k, k1, l, m, m2;
    double xr, xi;

    ip[0] = 0;
    l = n;
    m = 1;
    while ((m << 2) < l) {
        l >>= 1;
        for (j = 0; j <= m - 1; j++) {
            ip[m + j] = ip[j] + l;
        }
        m <<= 1;
    }
    if ((m << 2) > l) {
        for (i = 0; i <= n1 - 1; i++) {
            for (k = 1; k <= m - 1; k++) {
                for (j = 0; j <= k - 1; j++) {
                    j1 = (j << 1) + ip[k];
                    k1 = (k << 1) + ip[j];
                    xr = a[i][j1];
                    xi = a[i][j1 + 1];
                    a[i][j1] = a[i][k1];
                    a[i][j1 + 1] = a[i][k1 + 1];
                    a[i][k1] = xr;
                    a[i][k1 + 1] = xi;
                }
            }
        }
    } else {
        m2 = m << 1;
        for (i = 0; i <= n1 - 1; i++) {
            for (k = 1; k <= m - 1; k++) {
                for (j = 0; j <= k - 1; j++) {
                    j1 = (j << 1) + ip[k];
                    k1 = (k << 1) + ip[j];
                    xr = a[i][j1];
                    xi = a[i][j1 + 1];
                    a[i][j1] = a[i][k1];
                    a[i][j1 + 1] = a[i][k1 + 1];
                    a[i][k1] = xr;
                    a[i][k1 + 1] = xi;
                    j1 += m2;
                    k1 += m2;
                    xr = a[i][j1];
                    xi = a[i][j1 + 1];
                    a[i][j1] = a[i][k1];
                    a[i][j1 + 1] = a[i][k1 + 1];
                    a[i][k1] = xr;
                    a[i][k1 + 1] = xi;
                }
            }
        }
    }
}


void bitrv2row(int n, int n2, int *ip, double **a)
{
    int i, j, j1, k, k1, l, m;
    double xr, xi;

    ip[0] = 0;
    l = n;
    m = 1;
    while ((m << 1) < l) {
        l >>= 1;
        for (j = 0; j <= m - 1; j++) {
            ip[m + j] = ip[j] + l;
        }
        m <<= 1;
    }
    if ((m << 1) > l) {
        for (k = 1; k <= m - 1; k++) {
            for (j = 0; j <= k - 1; j++) {
                j1 = j + ip[k];
                k1 = k + ip[j];
                for (i = 0; i <= n2 - 2; i += 2) {
                    xr = a[j1][i];
                    xi = a[j1][i + 1];
                    a[j1][i] = a[k1][i];
                    a[j1][i + 1] = a[k1][i + 1];
                    a[k1][i] = xr;
                    a[k1][i + 1] = xi;
                }
            }
        }
    } else {
        for (k = 1; k <= m - 1; k++) {
            for (j = 0; j <= k - 1; j++) {
                j1 = j + ip[k];
                k1 = k + ip[j];
                for (i = 0; i <= n2 - 2; i += 2) {
                    xr = a[j1][i];
                    xi = a[j1][i + 1];
                    a[j1][i] = a[k1][i];
                    a[j1][i + 1] = a[k1][i + 1];
                    a[k1][i] = xr;
                    a[k1][i + 1] = xi;
                }
                j1 += m;
                k1 += m;
                for (i = 0; i <= n2 - 2; i += 2) {
                    xr = a[j1][i];
                    xi = a[j1][i + 1];
                    a[j1][i] = a[k1][i];
                    a[j1][i + 1] = a[k1][i + 1];
                    a[k1][i] = xr;
                    a[k1][i + 1] = xi;
                }
            }
        }
    }
}


void cftbcol(int n1, int n, double **a, double *w)
{
    int i, j, j1, j2, j3, k, k1, ks, l, m;
    double wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    for (i = 0; i <= n1 - 1; i++) {
        l = 2;
        while ((l << 1) < n) {
            m = l << 2;
            for (j = 0; j <= l - 2; j += 2) {
                j1 = j + l;
                j2 = j1 + l;
                j3 = j2 + l;
                x0r = a[i][j] + a[i][j1];
                x0i = a[i][j + 1] + a[i][j1 + 1];
                x1r = a[i][j] - a[i][j1];
                x1i = a[i][j + 1] - a[i][j1 + 1];
                x2r = a[i][j2] + a[i][j3];
                x2i = a[i][j2 + 1] + a[i][j3 + 1];
                x3r = a[i][j2] - a[i][j3];
                x3i = a[i][j2 + 1] - a[i][j3 + 1];
                a[i][j] = x0r + x2r;
                a[i][j + 1] = x0i + x2i;
                a[i][j2] = x0r - x2r;
                a[i][j2 + 1] = x0i - x2i;
                a[i][j1] = x1r - x3i;
                a[i][j1 + 1] = x1i + x3r;
                a[i][j3] = x1r + x3i;
                a[i][j3 + 1] = x1i - x3r;
            }
            if (m < n) {
                wk1r = w[2];
                for (j = m; j <= l + m - 2; j += 2) {
                    j1 = j + l;
                    j2 = j1 + l;
                    j3 = j2 + l;
                    x0r = a[i][j] + a[i][j1];
                    x0i = a[i][j + 1] + a[i][j1 + 1];
                    x1r = a[i][j] - a[i][j1];
                    x1i = a[i][j + 1] - a[i][j1 + 1];
                    x2r = a[i][j2] + a[i][j3];
                    x2i = a[i][j2 + 1] + a[i][j3 + 1];
                    x3r = a[i][j2] - a[i][j3];
                    x3i = a[i][j2 + 1] - a[i][j3 + 1];
                    a[i][j] = x0r + x2r;
                    a[i][j + 1] = x0i + x2i;
                    a[i][j2] = x2i - x0i;
                    a[i][j2 + 1] = x0r - x2r;
                    x0r = x1r - x3i;
                    x0i = x1i + x3r;
                    a[i][j1] = wk1r * (x0r - x0i);
                    a[i][j1 + 1] = wk1r * (x0r + x0i);
                    x0r = x3i + x1r;
                    x0i = x3r - x1i;
                    a[i][j3] = wk1r * (x0i - x0r);
                    a[i][j3 + 1] = wk1r * (x0i + x0r);
                }
                k1 = 1;
                ks = -1;
                for (k = (m << 1); k <= n - m; k += m) {
                    k1++;
                    ks = -ks;
                    wk1r = w[k1 << 1];
                    wk1i = w[(k1 << 1) + 1];
                    wk2r = ks * w[k1];
                    wk2i = w[k1 + ks];
                    wk3r = wk1r - 2 * wk2i * wk1i;
                    wk3i = 2 * wk2i * wk1r - wk1i;
                    for (j = k; j <= l + k - 2; j += 2) {
                        j1 = j + l;
                        j2 = j1 + l;
                        j3 = j2 + l;
                        x0r = a[i][j] + a[i][j1];
                        x0i = a[i][j + 1] + a[i][j1 + 1];
                        x1r = a[i][j] - a[i][j1];
                        x1i = a[i][j + 1] - a[i][j1 + 1];
                        x2r = a[i][j2] + a[i][j3];
                        x2i = a[i][j2 + 1] + a[i][j3 + 1];
                        x3r = a[i][j2] - a[i][j3];
                        x3i = a[i][j2 + 1] - a[i][j3 + 1];
                        a[i][j] = x0r + x2r;
                        a[i][j + 1] = x0i + x2i;
                        x0r -= x2r;
                        x0i -= x2i;
                        a[i][j2] = wk2r * x0r - wk2i * x0i;
                        a[i][j2 + 1] = wk2r * x0i + wk2i * x0r;
                        x0r = x1r - x3i;
                        x0i = x1i + x3r;
                        a[i][j1] = wk1r * x0r - wk1i * x0i;
                        a[i][j1 + 1] = wk1r * x0i + wk1i * x0r;
                        x0r = x1r + x3i;
                        x0i = x1i - x3r;
                        a[i][j3] = wk3r * x0r - wk3i * x0i;
                        a[i][j3 + 1] = wk3r * x0i + wk3i * x0r;
                    }
                }
            }
            l = m;
        }
        if (l < n) {
            for (j = 0; j <= l - 2; j += 2) {
                j1 = j + l;
                x0r = a[i][j] - a[i][j1];
                x0i = a[i][j + 1] - a[i][j1 + 1];
                a[i][j] += a[i][j1];
                a[i][j + 1] += a[i][j1 + 1];
                a[i][j1] = x0r;
                a[i][j1 + 1] = x0i;
            }
        }
    }
}


void cftbrow(int n, int n2, double **a, double *w)
{
    int i, j, j1, j2, j3, k, k1, ks, l, m;
    double wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    l = 1;
    while ((l << 1) < n) {
        m = l << 2;
        for (j = 0; j <= l - 1; j++) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            for (i = 0; i <= n2 - 2; i += 2) {
                x0r = a[j][i] + a[j1][i];
                x0i = a[j][i + 1] + a[j1][i + 1];
                x1r = a[j][i] - a[j1][i];
                x1i = a[j][i + 1] - a[j1][i + 1];
                x2r = a[j2][i] + a[j3][i];
                x2i = a[j2][i + 1] + a[j3][i + 1];
                x3r = a[j2][i] - a[j3][i];
                x3i = a[j2][i + 1] - a[j3][i + 1];
                a[j][i] = x0r + x2r;
                a[j][i + 1] = x0i + x2i;
                a[j2][i] = x0r - x2r;
                a[j2][i + 1] = x0i - x2i;
                a[j1][i] = x1r - x3i;
                a[j1][i + 1] = x1i + x3r;
                a[j3][i] = x1r + x3i;
                a[j3][i + 1] = x1i - x3r;
            }
        }
        if (m < n) {
            wk1r = w[2];
            for (j = m; j <= l + m - 1; j++) {
                j1 = j + l;
                j2 = j1 + l;
                j3 = j2 + l;
                for (i = 0; i <= n2 - 2; i += 2) {
                    x0r = a[j][i] + a[j1][i];
                    x0i = a[j][i + 1] + a[j1][i + 1];
                    x1r = a[j][i] - a[j1][i];
                    x1i = a[j][i + 1] - a[j1][i + 1];
                    x2r = a[j2][i] + a[j3][i];
                    x2i = a[j2][i + 1] + a[j3][i + 1];
                    x3r = a[j2][i] - a[j3][i];
                    x3i = a[j2][i + 1] - a[j3][i + 1];
                    a[j][i] = x0r + x2r;
                    a[j][i + 1] = x0i + x2i;
                    a[j2][i] = x2i - x0i;
                    a[j2][i + 1] = x0r - x2r;
                    x0r = x1r - x3i;
                    x0i = x1i + x3r;
                    a[j1][i] = wk1r * (x0r - x0i);
                    a[j1][i + 1] = wk1r * (x0r + x0i);
                    x0r = x3i + x1r;
                    x0i = x3r - x1i;
                    a[j3][i] = wk1r * (x0i - x0r);
                    a[j3][i + 1] = wk1r * (x0i + x0r);
                }
            }
            k1 = 1;
            ks = -1;
            for (k = (m << 1); k <= n - m; k += m) {
                k1++;
                ks = -ks;
                wk1r = w[k1 << 1];
                wk1i = w[(k1 << 1) + 1];
                wk2r = ks * w[k1];
                wk2i = w[k1 + ks];
                wk3r = wk1r - 2 * wk2i * wk1i;
                wk3i = 2 * wk2i * wk1r - wk1i;
                for (j = k; j <= l + k - 1; j++) {
                    j1 = j + l;
                    j2 = j1 + l;
                    j3 = j2 + l;
                    for (i = 0; i <= n2 - 2; i += 2) {
                        x0r = a[j][i] + a[j1][i];
                        x0i = a[j][i + 1] + a[j1][i + 1];
                        x1r = a[j][i] - a[j1][i];
                        x1i = a[j][i + 1] - a[j1][i + 1];
                        x2r = a[j2][i] + a[j3][i];
                        x2i = a[j2][i + 1] + a[j3][i + 1];
                        x3r = a[j2][i] - a[j3][i];
                        x3i = a[j2][i + 1] - a[j3][i + 1];
                        a[j][i] = x0r + x2r;
                        a[j][i + 1] = x0i + x2i;
                        x0r -= x2r;
                        x0i -= x2i;
                        a[j2][i] = wk2r * x0r - wk2i * x0i;
                        a[j2][i + 1] = wk2r * x0i + wk2i * x0r;
                        x0r = x1r - x3i;
                        x0i = x1i + x3r;
                        a[j1][i] = wk1r * x0r - wk1i * x0i;
                        a[j1][i + 1] = wk1r * x0i + wk1i * x0r;
                        x0r = x1r + x3i;
                        x0i = x1i - x3r;
                        a[j3][i] = wk3r * x0r - wk3i * x0i;
                        a[j3][i + 1] = wk3r * x0i + wk3i * x0r;
                    }
                }
            }
        }
        l = m;
    }
    if (l < n) {
        for (j = 0; j <= l - 1; j++) {
            j1 = j + l;
            for (i = 0; i <= n2 - 2; i += 2) {
                x0r = a[j][i] - a[j1][i];
                x0i = a[j][i + 1] - a[j1][i + 1];
                a[j][i] += a[j1][i];
                a[j][i + 1] += a[j1][i + 1];
                a[j1][i] = x0r;
                a[j1][i + 1] = x0i;
            }
        }
    }
}


void cftfcol(int n1, int n, double **a, double *w)
{
    int i, j, j1, j2, j3, k, k1, ks, l, m;
    double wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    for (i = 0; i <= n1 - 1; i++) {
        l = 2;
        while ((l << 1) < n) {
            m = l << 2;
            for (j = 0; j <= l - 2; j += 2) {
                j1 = j + l;
                j2 = j1 + l;
                j3 = j2 + l;
                x0r = a[i][j] + a[i][j1];
                x0i = a[i][j + 1] + a[i][j1 + 1];
                x1r = a[i][j] - a[i][j1];
                x1i = a[i][j + 1] - a[i][j1 + 1];
                x2r = a[i][j2] + a[i][j3];
                x2i = a[i][j2 + 1] + a[i][j3 + 1];
                x3r = a[i][j2] - a[i][j3];
                x3i = a[i][j2 + 1] - a[i][j3 + 1];
                a[i][j] = x0r + x2r;
                a[i][j + 1] = x0i + x2i;
                a[i][j2] = x0r - x2r;
                a[i][j2 + 1] = x0i - x2i;
                a[i][j1] = x1r + x3i;
                a[i][j1 + 1] = x1i - x3r;
                a[i][j3] = x1r - x3i;
                a[i][j3 + 1] = x1i + x3r;
            }
            if (m < n) {
                wk1r = w[2];
                for (j = m; j <= l + m - 2; j += 2) {
                    j1 = j + l;
                    j2 = j1 + l;
                    j3 = j2 + l;
                    x0r = a[i][j] + a[i][j1];
                    x0i = a[i][j + 1] + a[i][j1 + 1];
                    x1r = a[i][j] - a[i][j1];
                    x1i = a[i][j + 1] - a[i][j1 + 1];
                    x2r = a[i][j2] + a[i][j3];
                    x2i = a[i][j2 + 1] + a[i][j3 + 1];
                    x3r = a[i][j2] - a[i][j3];
                    x3i = a[i][j2 + 1] - a[i][j3 + 1];
                    a[i][j] = x0r + x2r;
                    a[i][j + 1] = x0i + x2i;
                    a[i][j2] = x0i - x2i;
                    a[i][j2 + 1] = x2r - x0r;
                    x0r = x1r + x3i;
                    x0i = x1i - x3r;
                    a[i][j1] = wk1r * (x0i + x0r);
                    a[i][j1 + 1] = wk1r * (x0i - x0r);
                    x0r = x3i - x1r;
                    x0i = x3r + x1i;
                    a[i][j3] = wk1r * (x0r + x0i);
                    a[i][j3 + 1] = wk1r * (x0r - x0i);
                }
                k1 = 1;
                ks = -1;
                for (k = (m << 1); k <= n - m; k += m) {
                    k1++;
                    ks = -ks;
                    wk1r = w[k1 << 1];
                    wk1i = w[(k1 << 1) + 1];
                    wk2r = ks * w[k1];
                    wk2i = w[k1 + ks];
                    wk3r = wk1r - 2 * wk2i * wk1i;
                    wk3i = 2 * wk2i * wk1r - wk1i;
                    for (j = k; j <= l + k - 2; j += 2) {
                        j1 = j + l;
                        j2 = j1 + l;
                        j3 = j2 + l;
                        x0r = a[i][j] + a[i][j1];
                        x0i = a[i][j + 1] + a[i][j1 + 1];
                        x1r = a[i][j] - a[i][j1];
                        x1i = a[i][j + 1] - a[i][j1 + 1];
                        x2r = a[i][j2] + a[i][j3];
                        x2i = a[i][j2 + 1] + a[i][j3 + 1];
                        x3r = a[i][j2] - a[i][j3];
                        x3i = a[i][j2 + 1] - a[i][j3 + 1];
                        a[i][j] = x0r + x2r;
                        a[i][j + 1] = x0i + x2i;
                        x0r -= x2r;
                        x0i -= x2i;
                        a[i][j2] = wk2r * x0r + wk2i * x0i;
                        a[i][j2 + 1] = wk2r * x0i - wk2i * x0r;
                        x0r = x1r + x3i;
                        x0i = x1i - x3r;
                        a[i][j1] = wk1r * x0r + wk1i * x0i;
                        a[i][j1 + 1] = wk1r * x0i - wk1i * x0r;
                        x0r = x1r - x3i;
                        x0i = x1i + x3r;
                        a[i][j3] = wk3r * x0r + wk3i * x0i;
                        a[i][j3 + 1] = wk3r * x0i - wk3i * x0r;
                    }
                }
            }
            l = m;
        }
        if (l < n) {
            for (j = 0; j <= l - 2; j += 2) {
                j1 = j + l;
                x0r = a[i][j] - a[i][j1];
                x0i = a[i][j + 1] - a[i][j1 + 1];
                a[i][j] += a[i][j1];
                a[i][j + 1] += a[i][j1 + 1];
                a[i][j1] = x0r;
                a[i][j1 + 1] = x0i;
            }
        }
    }
}


void cftfrow(int n, int n2, double **a, double *w)
{
    int i, j, j1, j2, j3, k, k1, ks, l, m;
    double wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    l = 1;
    while ((l << 1) < n) {
        m = l << 2;
        for (j = 0; j <= l - 1; j++) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            for (i = 0; i <= n2 - 2; i += 2) {
                x0r = a[j][i] + a[j1][i];
                x0i = a[j][i + 1] + a[j1][i + 1];
                x1r = a[j][i] - a[j1][i];
                x1i = a[j][i + 1] - a[j1][i + 1];
                x2r = a[j2][i] + a[j3][i];
                x2i = a[j2][i + 1] + a[j3][i + 1];
                x3r = a[j2][i] - a[j3][i];
                x3i = a[j2][i + 1] - a[j3][i + 1];
                a[j][i] = x0r + x2r;
                a[j][i + 1] = x0i + x2i;
                a[j2][i] = x0r - x2r;
                a[j2][i + 1] = x0i - x2i;
                a[j1][i] = x1r + x3i;
                a[j1][i + 1] = x1i - x3r;
                a[j3][i] = x1r - x3i;
                a[j3][i + 1] = x1i + x3r;
            }
        }
        if (m < n) {
            wk1r = w[2];
            for (j = m; j <= l + m - 1; j++) {
                j1 = j + l;
                j2 = j1 + l;
                j3 = j2 + l;
                for (i = 0; i <= n2 - 2; i += 2) {
                    x0r = a[j][i] + a[j1][i];
                    x0i = a[j][i + 1] + a[j1][i + 1];
                    x1r = a[j][i] - a[j1][i];
                    x1i = a[j][i + 1] - a[j1][i + 1];
                    x2r = a[j2][i] + a[j3][i];
                    x2i = a[j2][i + 1] + a[j3][i + 1];
                    x3r = a[j2][i] - a[j3][i];
                    x3i = a[j2][i + 1] - a[j3][i + 1];
                    a[j][i] = x0r + x2r;
                    a[j][i + 1] = x0i + x2i;
                    a[j2][i] = x0i - x2i;
                    a[j2][i + 1] = x2r - x0r;
                    x0r = x1r + x3i;
                    x0i = x1i - x3r;
                    a[j1][i] = wk1r * (x0i + x0r);
                    a[j1][i + 1] = wk1r * (x0i - x0r);
                    x0r = x3i - x1r;
                    x0i = x3r + x1i;
                    a[j3][i] = wk1r * (x0r + x0i);
                    a[j3][i + 1] = wk1r * (x0r - x0i);
                }
            }
            k1 = 1;
            ks = -1;
            for (k = (m << 1); k <= n - m; k += m) {
                k1++;
                ks = -ks;
                wk1r = w[k1 << 1];
                wk1i = w[(k1 << 1) + 1];
                wk2r = ks * w[k1];
                wk2i = w[k1 + ks];
                wk3r = wk1r - 2 * wk2i * wk1i;
                wk3i = 2 * wk2i * wk1r - wk1i;
                for (j = k; j <= l + k - 1; j++) {
                    j1 = j + l;
                    j2 = j1 + l;
                    j3 = j2 + l;
                    for (i = 0; i <= n2 - 2; i += 2) {
                        x0r = a[j][i] + a[j1][i];
                        x0i = a[j][i + 1] + a[j1][i + 1];
                        x1r = a[j][i] - a[j1][i];
                        x1i = a[j][i + 1] - a[j1][i + 1];
                        x2r = a[j2][i] + a[j3][i];
                        x2i = a[j2][i + 1] + a[j3][i + 1];
                        x3r = a[j2][i] - a[j3][i];
                        x3i = a[j2][i + 1] - a[j3][i + 1];
                        a[j][i] = x0r + x2r;
                        a[j][i + 1] = x0i + x2i;
                        x0r -= x2r;
                        x0i -= x2i;
                        a[j2][i] = wk2r * x0r + wk2i * x0i;
                        a[j2][i + 1] = wk2r * x0i - wk2i * x0r;
                        x0r = x1r + x3i;
                        x0i = x1i - x3r;
                        a[j1][i] = wk1r * x0r + wk1i * x0i;
                        a[j1][i + 1] = wk1r * x0i - wk1i * x0r;
                        x0r = x1r - x3i;
                        x0i = x1i + x3r;
                        a[j3][i] = wk3r * x0r + wk3i * x0i;
                        a[j3][i + 1] = wk3r * x0i - wk3i * x0r;
                    }
                }
            }
        }
        l = m;
    }
    if (l < n) {
        for (j = 0; j <= l - 1; j++) {
            j1 = j + l;
            for (i = 0; i <= n2 - 2; i += 2) {
                x0r = a[j][i] - a[j1][i];
                x0i = a[j][i + 1] - a[j1][i + 1];
                a[j][i] += a[j1][i];
                a[j][i + 1] += a[j1][i + 1];
                a[j1][i] = x0r;
                a[j1][i + 1] = x0i;
            }
        }
    }
}

double **alloc_2d_double(int n1, int n2)
{
    double **dd, *d;
    int j;

    dd = (double **) malloc(sizeof(double *) * n1);
    alloc_error_check(dd);
    d = (double *) malloc(sizeof(double) * n1 * n2);
    alloc_error_check(d);
    dd[0] = d;
    for (j = 1; j < n1; j++) {
        dd[j] = dd[j - 1] + n2;
    }
    return dd;
}

void free_2d_double(double **dd)
{
    free(dd[0]);
    free(dd);
}

int *alloc_1d_int(int n1)
{
    int *i;

    i = (int *) malloc(sizeof(int) * n1);
    alloc_error_check(i);
    return i;
}

void free_1d_int(int *i)
{
    free(i);
}

double *alloc_1d_double(int n1)
{
    double *d;

    d = (double *) malloc(sizeof(double) * n1);
    alloc_error_check(d);
    return d;
}

void free_1d_double(double *d)
{
    free(d);
}
