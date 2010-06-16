#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QProcess>
#include <QInputDialog>


class MainWindow : public QInputDialog
{
  Q_OBJECT
  public:
  MainWindow(QWidget* parent=0);
 
public slots:
  void postProcess();
  void check(const QString&);
  void vcheck();
  // void help(const QString&);
  void vmerge(); //vogmerge button clicked
  void fvmAdapt(); //FVMAdapt
  void import();
  void showQuality(QString, QProcess::ExitStatus, QString);
  void processItem(const QString&);
  void generateVar();
signals:
 
private:
  bool waitForQualityFile;
};

#endif

