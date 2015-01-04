#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "DynamicArray2D.h"

namespace Ui {
class MainWindow;
}

class BmpLabel;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    DynArray2D<double> kinoformArray;  // キノフォーム
    DynArray2D<double> intensityArray; // intensity distribution

    int kinoformAreaSize;                //生成するキノフォームのサイズ
    int focusNumber;                     //焦点数
    int convertCount;

    // fourier transform
    void cdft2d(int, int, int, double **, int *, double *);
    void laserIntensity(double, double);
    //E:物体面の電場→高速フーリエ変換→ホログラム面の電場
    //W:三角関数テーブル
    int *ip;
    double **E, **E0, *W;

private slots:
    void on_showIntensityButton_clicked();
    void on_showKinoformButton_clicked();
    void on_convertButton_clicked();
    void on_saveButton_clicked();

    void getUIWidgetValue();
    void updateUIWidget();

private:
    Ui::MainWindow *ui;
    BmpLabel *confirmBmpLabel_;

    double waveLength_;                  //波長
    double focalLength_;                 //焦点距離
    double beamWaist_;                   //ビームウェスト
    int	tc_;                             //トポロジカルチャージ
    double initialPhase_l_;              //光渦の初期位相
    double r_;                           //焦点面における焦点移動量
    double phi_;                         //焦点面における焦点移動方向
    double initialPhase_r_;              //焦点面における焦点移動の初期位相
    double z_;                           //光軸上の焦点移動量
    double initialPhase_z_;              //光軸上の焦点移動における初期位相

private:
    void init();
    void initKinoform();
    void outputIntensity();
    void outputKinoform();

    // kernel function
    void appliedLGBeam();
    void appliedFlesnel();
    void appliedFocalShift();
};

#endif // MAINWINDOW_H
