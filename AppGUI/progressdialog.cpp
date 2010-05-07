

#include <QProcess>
#include <QTimer>
#include <QtDebug>
#include <QLabel>
#include <QVBoxLayout>
#include <QPushButton>
#include <QTextEdit>
#include <QDir>
#include "progressdialog.h"








ProgressDialog::ProgressDialog(QString command, bool autoclose, QWidget* parent):QWidget(parent), autoClose(autoclose){

  
  QWidget::setAttribute(Qt::WA_DeleteOnClose);


  commands = command.split('\n');
  for(int i = 0; i < commands.size(); i++){
    args<< commands[i].section(' ', 1, -1).split(' ');
    programs << commands[i].section(' ', 0, 0);
    proc << 0;

  }
  // qDebug() << programs;
  //qDebug() <<args;
  
  nameLabel = new QLabel("      ");
  timeLabel = new QLabel("0s");
  text = new QTextEdit;
  abortButton= new QPushButton("&Abort");
  
 
  QHBoxLayout* labelLayout = new QHBoxLayout;
  labelLayout->addWidget(nameLabel);
  labelLayout->addWidget(timeLabel);
  labelLayout->addWidget(abortButton);
  QVBoxLayout* mainLayout = new QVBoxLayout;
  mainLayout->addLayout(labelLayout);
  mainLayout->addWidget(text);
  setLayout(mainLayout);  
  
  
  if(programs.size() > 0)startProgress();
 
}




bool ProgressDialog::startProgress()
{
  proc[0] = new QProcess;
  t = new QTimer(this);
  currentP = 0;
  connect(proc[0], SIGNAL(finished(int, QProcess::ExitStatus)),
          this, SLOT(next(int, QProcess::ExitStatus)));
  

  connect(proc[0], SIGNAL(readyReadStandardError()), this, SLOT(setText()));
  connect(proc[0], SIGNAL(readyReadStandardOutput()), this, SLOT(setText()));
  
  connect(t, SIGNAL(timeout()), this, SLOT(update()));
  connect(abortButton, SIGNAL(clicked()), proc[0], SLOT(kill())); 
  
  proc[0]->setWorkingDirectory(QDir::currentPath());
  proc[0]->setProcessChannelMode(QProcess::MergedChannels);  
  if (proc[0]->state() != QProcess::Running) {
    proc[0]->start(programs[0], args[0]);
    if (!proc[0]->waitForStarted()) {
      nameLabel->setText("can not launch " + commands[0]);
    }
  }
  //start timeing
  {
    t->start(1000);
    time.start();
    nameLabel->setText(commands[0]+ " running");
  }
  
  return true;
}

void ProgressDialog::updateState(QProcess::ProcessState status){
  switch(status){
  case QProcess::Running:
    t->start(1000);
    time.start();
    nameLabel->setText(commands[currentP]+ " running");
    break;
  case QProcess::NotRunning:
    nameLabel->setText(commands[currentP]+" finished");
    t->stop();
    
    break;
  case QProcess::Starting:
    nameLabel->setText("starting " );
    break;
  }
}

void ProgressDialog::update(){
  timeLabel->setText(QString("%1s").arg(time.elapsed()/1000));
}
void ProgressDialog::setText(){
  QByteArray output = proc[currentP]->readAllStandardOutput();
  text->insertPlainText(QString(output));
}

void ProgressDialog::next(int exitCode, QProcess::ExitStatus exitStatus){
  t->stop();
  proc[currentP]->deleteLater();
  
  
  if(exitStatus==QProcess::CrashExit || exitCode != 0){
   
    nameLabel->setText(commands[currentP]+ " crashed") ;
    emit progressFinished(commands[currentP], QProcess::CrashExit);
    //  if(autoClose)close(); 
    return; 
  }else if(currentP == (programs.size()-1)){
   
    nameLabel->setText(commands[currentP]+ " exit normally");
    emit progressFinished(commands[currentP], QProcess::NormalExit);
    if(autoClose)close();
    return;
  }
  
  
  currentP++;
  if(currentP == programs.size()) return;

  proc[currentP] = new QProcess;
  connect(proc[currentP], SIGNAL(finished(int, QProcess::ExitStatus)),
          this, SLOT(next(int, QProcess::ExitStatus)));
  
  // connect(proc, SIGNAL(stateChanged(QProcess::ProcessState)),
  //           this, SLOT(updateState(QProcess::ProcessState)));
  connect(proc[currentP], SIGNAL(readyReadStandardError()), this, SLOT(setText()));
  connect(proc[currentP], SIGNAL(readyReadStandardOutput()), this, SLOT(setText()));
  
  connect(t, SIGNAL(timeout()), this, SLOT(update()));
  connect(abortButton, SIGNAL(clicked()), proc[currentP], SLOT(kill()));
  
  proc[currentP]->setWorkingDirectory("/simcenter/data1/qxue/grids/");
  proc[currentP]->setProcessChannelMode(QProcess::MergedChannels);
   
  if (proc[currentP]->state() != QProcess::Running) {
    proc[currentP]->start(programs[currentP], args[currentP]);
    if (!proc[currentP]->waitForStarted()) {
      nameLabel->setText("can not launch " + commands[currentP]);
    }else{
      nameLabel->setText(commands[currentP]+" running");
    }
  }
  t->start(1000);
  time.start(); 
  
 }

