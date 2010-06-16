#ifndef PROGRESSDIALOG_H
#define PROGRESSDIALOG_H

#include <QProgressDialog>
#include <QProcess>
#include <QStringList>
#include <QList>
#include <QString>
#include <QWidget>
#include <QTime>

class QTimer;

class QWidget;
class QLabel;
class QPushButton;
class QTextEdit;
class ProgressDialog: public QWidget
{
  Q_OBJECT
  
  public:
  ProgressDialog(QString command,QString directory, bool autoclose = true, QWidget *parent = 0);
  //  ~ProgressDialog();
 
signals:
  void progressFinished(QString, QProcess::ExitStatus, QString);
public slots:
  void updateState(QProcess::ProcessState);
  void next(int, QProcess::ExitStatus);
  void update();
  void setText();
private:
   bool startProgress();

  QString workingDirectory;
  QStringList commands;
  QLabel* nameLabel;
  QLabel* timeLabel;
  QTextEdit* text;
  QPushButton* abortButton;
  QList<QProcess*> proc;
  QStringList  programs;
  QList<QStringList> args;
  QTimer* t;
  int currentP;
  QTime time;
  bool autoClose;
};
#endif  
