#include <QtGui>
#include <QDomElement>
#include "getfile.h"

GetFileWindow::GetFileWindow(QString exp, QString& fileSelected, QWidget *parent) : QDialog(parent)
 {
     browseButton = createButton(tr("&Browse..."), SLOT(browse()));
     findButton = createButton(tr("&Find"), SLOT(find()));

     fileComboBox = createComboBox(exp);
     //     textComboBox = createComboBox();
     directoryComboBox = createComboBox(QDir::currentPath()+"/");

     fileLabel = new QLabel(tr("Named:"));
     // textLabel = new QLabel(tr("Containing text:"));
     directoryLabel = new QLabel(tr("In directory:"));
     filesFoundLabel = new QLabel;

     createFilesTable();
     fileNameLabel = new QLabel;
    

       QHBoxLayout *buttonsLayout = new QHBoxLayout;
      buttonsLayout->addStretch();
      buttonsLayout->addWidget(findButton);
      //  buttonsLayout->setSpacing(40);
      buttonsLayout->addWidget(fileNameLabel);
      buttonsLayout->setAlignment(findButton, Qt::AlignLeft);
       buttonsLayout->setAlignment(fileNameLabel, Qt::AlignRight);
     QGridLayout *mainLayout = new QGridLayout;
     
   
     mainLayout->addWidget(fileLabel, 0, 0);
     
     mainLayout->addWidget(fileComboBox, 0, 1, 1, 2);
     // mainLayout->addWidget(textLabel, 1, 0);
     // mainLayout->addWidget(textComboBox, 1, 1, 1, 2);
     mainLayout->addWidget(directoryLabel, 1, 0);
     mainLayout->addWidget(directoryComboBox, 1, 1);
     mainLayout->addWidget(browseButton, 1, 2);
     mainLayout->addWidget(filesTable, 2, 0, 1, 3);
     mainLayout->addWidget(filesFoundLabel, 3, 0);
     mainLayout->addLayout(buttonsLayout, 4, 0);
     //mainLayout->addWidget(fileNameLabel, 4, 1, 1, 20);
     setLayout(mainLayout);
     fileNameLabel->hide();
     selectedFileName = fileSelected ;
     setWindowTitle(tr("Find Files"));
     findButton->setMaximumWidth(40);
     resize(700, 300);
 }
void GetFileWindow::addDirectory(QString dir){
  QStringList dir_list = dir.split(",", QString::SkipEmptyParts);
  directoryComboBox->addItems(dir_list);
  
}
 void GetFileWindow::browse()
 {
     QString directory = QFileDialog::getExistingDirectory(this,
                                tr("Find Files"), QDir::currentPath());
     if (!directory.isEmpty()) {
       int index = directoryComboBox->findText(directory);
       if(index != -1)directoryComboBox->setCurrentIndex(index);
       else{
         directoryComboBox->addItem(directory);
         directoryComboBox->setCurrentIndex(directoryComboBox->count()-1);
       }
     }
 }

 void GetFileWindow::find()
 {
     filesTable->setRowCount(0);

     QString fileName = fileComboBox->currentText();
     // QString text = textComboBox->currentText();
     QString path = directoryComboBox->currentText();

     QDir directory = QDir(path);
     QStringList files;
     if (fileName.isEmpty())
         fileName = "*";
     files = directory.entryList(QStringList(fileName),
                                 QDir::Files | QDir::NoSymLinks);

     // if (!text.isEmpty())
     //  files = findFiles(directory, files, text);
     showFiles(directory, files);
 }

 QStringList GetFileWindow::findFiles(const QDir &directory, const QStringList &files,
                               const QString &text)
 {
     QProgressDialog progressDialog(this);
     progressDialog.setCancelButtonText(tr("&Cancel"));
     progressDialog.setRange(0, files.size());
     progressDialog.setWindowTitle(tr("Find Files"));

     QStringList foundFiles;

     for (int i = 0; i < files.size(); ++i) {
         progressDialog.setValue(i);
         progressDialog.setLabelText(tr("Searching file number %1 of %2...")
                                     .arg(i).arg(files.size()));
         qApp->processEvents();

         if (progressDialog.wasCanceled())
             break;

         QFile file(directory.absoluteFilePath(files[i]));

         if (file.open(QIODevice::ReadOnly)) {
             QString line;
             QTextStream in(&file);
             while (!in.atEnd()) {
                 if (progressDialog.wasCanceled())
                     break;
                 line = in.readLine();
                 if (line.contains(text)) {
                     foundFiles << files[i];
                     break;
                 }
             }
         }
     }
     return foundFiles;
 }

 void GetFileWindow::showFiles(const QDir &directory, const QStringList &files)
 {
     for (int i = 0; i < files.size(); ++i) {
         QFile file(directory.absoluteFilePath(files[i]));
         qint64 size = QFileInfo(file).size();

         QTableWidgetItem *fileNameItem = new QTableWidgetItem(files[i]);
         fileNameItem->setFlags(fileNameItem->flags() ^ Qt::ItemIsEditable);
         QTableWidgetItem *sizeItem = new QTableWidgetItem(tr("%1 KB")
                                              .arg(int((size + 1023) / 1024)));
         sizeItem->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
         sizeItem->setFlags(sizeItem->flags() ^ Qt::ItemIsEditable);

         int row = filesTable->rowCount();
         filesTable->insertRow(row);
         filesTable->setItem(row, 0, fileNameItem);
         filesTable->setItem(row, 1, sizeItem);
     }
     filesFoundLabel->setText(tr("%1 file(s) found").arg(files.size()) +
                              (" ( Click on a file to select it)"));
 }

 QPushButton *GetFileWindow::createButton(const QString &text, const char *member)
 {
     QPushButton *button = new QPushButton(text);
     connect(button, SIGNAL(clicked()), this, member);
     return button;
 }

 QComboBox *GetFileWindow::createComboBox(const QString &text)
 {
     QComboBox *comboBox = new QComboBox;
     comboBox->setEditable(true);
     comboBox->addItem(text);
     comboBox->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
     return comboBox;
 }

 void GetFileWindow::createFilesTable()
 {
     filesTable = new QTableWidget(0, 2);
     filesTable->setSelectionBehavior(QAbstractItemView::SelectRows);

     QStringList labels;
     labels << tr("File Name") << tr("Size");
     filesTable->setHorizontalHeaderLabels(labels);
     filesTable->horizontalHeader()->setResizeMode(0, QHeaderView::Stretch);
     filesTable->verticalHeader()->hide();
     filesTable->setShowGrid(false);

     connect(filesTable, SIGNAL(cellClicked(int, int)),
             this, SLOT(openFileOfItem(int, int)));
 }


 void GetFileWindow::openFileOfItem(int row, int /* column */)
 {
     QTableWidgetItem *item = filesTable->item(row, 0);
     selectedFileName =  directoryComboBox->currentText() + item->text();
    
     fileNameLabel->setText("selected file: " + selectedFileName);
     fileNameLabel->show();
     emit fileNameSelected(selectedFileName);
     //  QDesktopServices::openUrl(item->text());
 }

FindFileWindow::FindFileWindow(const QDomElement& my_elem, QString& fileSelected, QWidget *parent)
  : GetFileWindow(my_elem.attribute("exp"), fileSelected, parent){
  myelem=my_elem;
  if(myelem.hasAttribute("current"))emit fileNameSelected(myelem.attribute("current"));
  connect(this, SIGNAL(fileNameSelected(QString)), this, SLOT(updateElem(QString)));
}
void FindFileWindow::updateElem(QString s){
  myelem.setAttribute("current",s);
}
