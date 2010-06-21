#ifndef PAGES_H
#define PAGES_H
#include <QWidget>
#include <QDomDocument>
#include <QString>
#include <QObject>
#include <QDomElement>
#include <QLineEdit>
#include <QGroupBox>
#include <QStringList>
#include <QDoubleValidator>
#include <QSignalMapper>
#include <QList>
#include <QComboBox>
#include <QPointer>
#include <iostream>
#include <QString>
#include <QItemDelegate>
#include <QDoubleSpinBox>
#include "grid.h"
class QPushButton;
class QButtonGroup;
class QLabel;
class QListWidget;
class QStackedWidget;
class QListWidgetItem;
class QCheckBox;


bool angleBetween(positions3d v1,positions3d v2, double& heading, double& attitude, double& bank);


bool conditionIsSatisfied(const QDomElement& theroot, const QString& condition);






class DoubleEdit : public QLineEdit{
  
  Q_OBJECT
  
  public:
  DoubleEdit(QWidget *parent = 0);
  DoubleEdit(double d, QWidget *parent = 0);
  void setBottom(double);
  void setTop(double);
  void setRange(double minn, double maxx);
  double value();
public slots:
  void mapValue(int);
  void setValue(double);
private slots:
  void changeValue(const QString&);
signals:
  void valueChanged(double);
private:
  QDoubleValidator* validator;
};



class DoubleSpinBox:public QDoubleSpinBox{
  Q_OBJECT
  public:
  DoubleSpinBox(QWidget* parent = 0);
signals:
  void paraChanged();//this signal is emitted when the step value or decimal value changed
protected:
  void keyPressEvent(QKeyEvent *event);
};



class LabeledDoubleSpBox : public QWidget{
  
  Q_OBJECT
  
  public:
  
  LabeledDoubleSpBox(const QString& title,  QWidget *parent = 0);
  void setRange(double minn, double maxx);
  double value();
  inline double singleStep(){return edit->singleStep();}
  inline int decimals(){return edit->decimals();}
  inline void setSingleStep(double d){edit->setSingleStep(d);}
  inline void setDecimals(int n){edit->setDecimals(n);} 
public slots:
  void setValue(double);
  void display();
  void undisplay();
private slots:
  
signals:
  void valueChanged(double);
  void paraChanged();
protected:
  DoubleSpinBox* edit;
  QLabel* valueLabel;
};


  
class VectSpBox:public QGroupBox{
  Q_OBJECT
  
  public:
    
  VectSpBox( const QString& title, QWidget *parent = 0);
  void setRange(double d1, double d1);
  void setXRange(double minn, double maxx);
  void setYRange(double minn, double maxx);
  void setZRange(double minn, double maxx);
  positions3d value();
public slots:
  void setValue(const positions3d& p);
private slots:
  void setInfo();
  void setPara();
signals:
  void valueChanged(const positions3d&);
private:
  LabeledDoubleSpBox* xedit;
  LabeledDoubleSpBox* yedit;
  LabeledDoubleSpBox* zedit;
 
};


class IntEdit : public QLineEdit{
  
  Q_OBJECT
  
  public:
  IntEdit(QWidget *parent = 0);
  IntEdit(int d, QWidget *parent = 0);
  void setBottom(int);
  void setTop(int);
  void setRange(int minn, int maxx);
 
  int value();
public slots:
  void setValue(int);
private slots:
  void changeValue(const QString&);
signals:
  void valueChanged(int);
private:
  QIntValidator* validator;
};


class GeneralGroup : public QGroupBox
{
  Q_OBJECT
  
  public:
  GeneralGroup(QDomElement& elem,QWidget* parent = 0);
  QString currentText();
public slots:
  void updateComponents(); //for dvector
  void changeState();
  void updateChecked();
  void updateShowStatus(const bool&);
  void gotoUnfinishedPage(){};
signals:
  void updateStatus(const QString&);
  void updateStatusTip(int);
  void stateChanged();
  void componentsChanged();
  void textChanged(const QString &);
  void showStatus(const bool &);
protected:
  QDomElement myelem;

};




class VarGBox : public GeneralGroup{
  
  Q_OBJECT
  
  public:
  VarGBox(QDomElement& elem, QWidget *parent = 0);
  QString currentText();
signals:
 
public slots:
  void updateComponents(); //for dvector
  void changeState();
private slots:
  void updateCurrentX(const QString&);
  void updateCurrentY(const QString&);
  void updateCurrentZ(const QString&);
  void updateCurrentUnit(const QString&);
  void updateUnitList(const QString&);
  void updateLabels(int );
  void updateCurrent(const QString&);
  void updateShowStatus(const bool&);
  void updateSelection(int); //for flags
  void setDefault(bool satisfied);//for default condition, when the condition is satisfied

private:
 /*!
   for vectors and dvectors 
 */ 
  QList<QLabel*> labels;
  QList<QWidget*> mfs;
 
};

class AllGroup : public GeneralGroup{
  
  Q_OBJECT
  
  public:
  AllGroup(  QDomElement& elem,  bool isVertical = true, QWidget *parent = 0);
  QString currentText();
  
signals:
public slots:
  void childrenClicked(int i);//for 2of3 element
  void updateCurrentText();
 
private:
  //for 2of3 element
  QList<QWidget*> objs;
  QSignalMapper *signalMapper;
};


  

class StackGroup : public GeneralGroup{
  
  Q_OBJECT
  
  public:
  StackGroup(  QDomElement& elem, QWidget *parent = 0);
  QString currentText();
  void clearText();
signals:
  
public slots:
  void updateCurrentText();
  void changePage(int);
private:
  
  QComboBox *typesWidget;
  QStackedWidget *pagesWidget;
  QStringList whatsThisList;
  QStringList toolTipList;

};






//class for a panel
class VarPanel : public GeneralGroup{
  Q_OBJECT
  public:
  VarPanel(QDomElement &elem,
          QWidget *parent=0);
  ~VarPanel();
  QString currentText();
  void clearText();
  
public slots:
  void advancedButtonClicked(int);
  void updateCurrentText();
  void updateShowStatus(const bool&);
signals:

  void updateStatusTip(int);
  void updateStatus( QString);

protected:
  QList<QWidget*> myAdvancedGroup;
  QButtonGroup* buttonGroup;
  QLabel* currentLabel;
  bool showStatus;
};



//for physicswindow and solverwindow


class AllVBWindow : public GeneralGroup{
  
  Q_OBJECT
   
  public:
  AllVBWindow(  QDomElement& elem, QWidget *parent = 0);
  
public slots:
  bool save();
private:

};

class ChoiceGroup : public GeneralGroup{
  Q_OBJECT
  public:
  ChoiceGroup( QDomElement &elem,
               QWidget *parent=0);
  ~ChoiceGroup();
  QString currentText();
public slots:
  void editButtonPressed();
  void changeState();
  void update(int);
  void updateCurrentText();
  void updateShowStatus(const bool&);
signals:
private:
  
  QPushButton* editButton;
  QWidget* editGroup;
  QComboBox* comboBox;
  QStringList editItems;
  QStringList whatsThisList;
  QStringList toolTipList;
  bool showStatus;
 
};


class StateStackGroup : public GeneralGroup{
  
  Q_OBJECT
  
  public:
  StateStackGroup(  QDomElement& elem, 
                QWidget *parent = 0);
  QString currentText();
  
signals:
   /*!
     This signal is emited when my state changed. To tell \a parent update its state
   */ 
  void stateChanged(QString);

  
public slots:
  void changePage(int i);//when buttonGroup clicked
  QString myState();
  void add();
  void updateCurrentText();
  

private:
  QButtonGroup *buttonGroup;
  QStackedWidget *pagesWidget;
  
};






class VarPage : public GeneralGroup{
  Q_OBJECT
  public:
  VarPage(QDomElement &elem,
       QWidget *parent=0);
   
private slots:
  void updateCurrentText();
signals:
  void textChanged(const QString &);
protected:
};

// Class that delegates for the boundary visibility model
class showDelegate: public QItemDelegate
{
  Q_OBJECT
  public:
  QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem &option,
			const QModelIndex &index) const;
  showDelegate(QObject* parent=0):QItemDelegate(parent){};
  
};

class colorDelegate: public QItemDelegate
{
  
  Q_OBJECT
  public:
  bool  editorEvent(QEvent *event, QAbstractItemModel *model,
                    const QStyleOptionViewItem &option,
                    const QModelIndex &index);
  colorDelegate(QObject* parent=0):QItemDelegate(parent){};
  
};





#endif
/*! \mainpage  chemdemo documentation

  \ref page_xml
  
    
  
  \page page_xml how to use xml file to modify chemdemo GUI
  
  \section whyxml  Why xml file is used
 
  Chem is a solver for fluid dynamics simulation, and Chemdemo is a GUI
  designed for the beginners to generate Chem input variable file conveniently. Like
  any other software that is within its life cycle,  Chem is a tool that
  still under evolution. In order for Chemdemo to grow with Chem,
  xml files is used so that the seasoned users can modify GUI by
  changing xml files.
 

  Not like qt .ui file, the xml files here are not focused on the geometry
  design, even though the geometry layout may be mentioned.  It's used to
  describe the properties and the grouping of the input variables.
  
  
  In the following part of this document, the semantics/syntax of these xml
  files will be described in a bottom-up way.  
  

  \section xml_variable How to define an input variable in xml file
  
  
  Since the input .var file of Chem is the values of a set of variables,
  the properties of these variables need to be defined in xml files. To make the xml more
  readable and easier to modify, the way the properties are defined
  should be consistent and instinctive.
  
  In xml file, the commonly used variables are defined in Dom element \a variables, others are defined in the context where they are used. If in the place where a variable is used, the attribute \a element is not specified, the program will go to \a variables to search for the variable,  then merge the attributes into the element where the varaible is used.
  
  Each variable is represented by one DOM element in  xml file, and
  each property of the variable is described by an attribute of the DOM element.  


  

 
  Here are some attributes that can be used to define properties of variables in xml file:
  
 
-#  \b element: element is used to describe the data type of the variable.
       The values of attribute  `element' can be: float,
       int, string, selection, vector and dvector. This attribute need to be
       specified for all variables.\n\n
       For `float' and `int' data types,
       attributes `bottom',  `top' and `range' can be used to describe the
       minimum value, the maximum value and the range of this variable.\n\n
       `selection' data type is used for flags. One or more flags can be
       selected, and all the choices are defined as  child elements of the
       variable. Here is the example for `selection' data type:\n\n
       < flags   element="selection" >\n
       <normal/>\n
       <equilibrium/>\n
       </flags>\n\n
       `vector' data type is used for 3d vector such as velocity. In GUI,
       the user can select if the vector is represented in Cartesian
       coordinates system or Polar coordinates system. But these details
       are hard coded and  are hidden from xml files\n\n
-# \b default:  the default value of the variable. Used when attribute 'defaultCondition' is defined. \n\n
-# \b name: the attribute `name' is the variable name that appears  in .vars
     file. For example, pressure in .vars file is represented as `p', `p0'
     or `PWall'. And `p', `p0', `pwall' are the possible values of
     attribute `name' for variable pressure. Not all variables have names,
     flags only have  values.  \n\n
-# \b prefix: in .vars file, the symbol between variable name and its
     value. The default prefix is `='. Sometimes it's ` : '. \n\n 
-#  \b title: title is the name or label that appear on the interface. Most
     of the time, title doesn't need to be specified. If title is not
     specified, the tag name of the variable will be used as the title.\n\n 
-#  \b bottom, \b top, \b range: the attribute for int or float data type. As
     introduced in attribute `element'\n\n
-#  \b unit: the attribute of unit is a string list separated by comma. If
     this attribute is defined, a combo box will appear on the interface,
     display all the possible units of the variable. For example, the
     attribute  `unit' of pressure is:\n
     unit="atm,Pa,psi,lbf/in/in,lbf/ft/ft,slug/ft/ft"\n
     If a variable has a unit, the attribute `unit' need to be specified. And the first unit
     in the list is the default unit. The user can add more units by
     editing from combo box.\n\n
-#  \b condition: some variables only need to be specified if  a condition is
     satisfied. These conditions are the physics state that the system is
     in, such as whether the system is steady or not, whether the system has
     single component or multiple components.\n\n
     To specify  the conditions of a variable, the user need look up the
     \a physicsPage in xml file, and use the equivalence between
     tag names
     to define the conditions. For example, if a variable is only needed in
     a multiple component system, the attribute `condition' can be defined
     as: \n
     condition=''Material==multi''\n\n
-#  \b defaultCondition: There are two kinds of conditions. One is that if the condition is
     satisfied, the variable will exist; otherwise, the variable will not
     exist. We use attribute `condition' to define it. The other
     kind of condition is that if the condition is not satisfied, it can be
     any value within its range, otherwise, it can only be a default
     value. We use attribute `defaultCondition' and `default' together to
     specify it.  
   
 
 

The following attribute are used to represent the status of the
variables. Its value can be modified by the program, and usually the
user specify it only when it's necessary.      

-# \b current: the current value.\n\n
-# \b conditionSatisfied:  specified with `condition' or
    `defaultCondition'. It represents if the condition is satisfied at
    this moment.\n\n
-# \b checked: specified only if the variable is  optional. The values of
    attribute can be `1' or `0'. `1' means the the variable is selected,
    `0' means it's not.

Some attributes are used to provide help information, these attributes
also applies to other elements that represented by a QWidget.


-# \b toolTip: The toolTip is a short piece of text reminding the user of the
    widget's function.It is drawn immediately below the given position
    in a distinctive black-on-yellow color combination. When the mouse
    is on the widget, the text will be drawn.\n\n 
-# \b whatsThis:  help texts are typically longer and more detailed than
    tooltips, but generally provide less information than that supplied by
    separate help windows.\n\n
    QWhatsThis provides a single window with an explanatory text that pops
    up when the user asks "What's This?". The default way for users to ask
    the question is to move the focus to the relevant widget and press
    Shift+F1. The help text appears immediately; it goes away as soon as
    the user does something else.\n\n
-# \b statusTip:  a short piece of text displayed in statusbar temporary
    when the mouse is in the area of the widget.\n\n 

\sa \ref VarGBox
\section variableGroup Grouping of variables

 The Grouping of variable is putting a set of related variables into one
  group box. The attribute `element' is used to define the relation
  between variables, its values  can be: all, 2of3, choice or
  stack. The child of a group can be a group or a variable.

  
-# \b all: all the children need to be specified.\n\n
-# \b stack: only one of the children need to be specified. On the
    interface, There is an exclusive radio button group, if the user checks one
    of them, the interface of the selected variable will be brought up. \n\n
-# \b 2of3: the group has three children, only two of them need to be
    specified.\n\n  
-# \b choice:  actually, only one variable is defined.   And it's a selection of one value from many values. All
    choices are listed as child elements. For example,\n\n
    <transport_model  element="choice" >\n
    <sutherland  editable="Sland" />\n
    <powerLaw  editable="powerLawParam" />\n
    <chemkin   />\n
    <const_viscosity  editable="const_viscosity" />\n
    </transport_model>\n\n
    However, `choice' type allows its children to use attribute
   `editable' to define another compound variable. So we consider
   variable with element `choice' as a group. In the former example,
   compound variable Sland and const_viscosity can be found in the
   element \a editable.    

\section xml_panel Panel

A panel is used to specify a set of variables or a compound
variable. The compound variable is a variable that consists of a set
of variables, such as  initialConditions, inflow, etc. A
panel can have groups and variables as its children. It's different
from a group with element=`all' in three ways:


-  Advanced features: It may have child element \a advanced. We can consider \a advanced' as a
   special kind of group for hiding the
   advanced featured from the beginners. On GUI, when the button \a advanced is pressed,  the interface of
   the advanced features will pop up as a
   new window.\n\n 
-  Feedback features: There is a label at the
   bottom of the panel. Whenever the user changes any input from the
   panel, the current text displayed on the label will change
   instantly. This feature can be disabled to save space.\n\n
-  Output as a compound variable:
   when a group is written into .vars files, all its childrens are
   written out independently. A panel might be written out as one compound
   variable if the attribute \a prefix is defined as a nonempty string.
   For example, Sland:<a1=1.2e-04, a2=...>, or inflow(p=1atm, T=300K ...)
  






The following attributes define the format of the text of panels: 


-# \b prefix: the symbols between the name of the variable and its value. The
    default value is empty string, which means the panel is used to
    specify a set of independent variables instead of a compound
    variable. If \prefix is not empty, a compound variable is defined.
    For example, for initial condition, prefix=`:< '.  And for boundary condition,  prefix=`('.\n\n
-# \b postfix:  the symbols that are appended after the text of last
    child. Default value: empty string.\n\n
-# \b infix:  the symbols that are inserted between the text of child
    elements. The default value is newline character.
 
  
  
The following attributes define the layout of the panel:


- \b numColumn: the panels use  QGridLayout. By default, numColumn=1,
    i.e., all the children are aligned vertically. If numColumn >=
    number-of-children, the panel is aligned horizontally. otherwise,  1<
    numColumn< number-of-children, the panel has a grid layout.  


More layout features will be added in the future.
The other attributes of a panel include:


-  \b feedback:  By default, the instant feedback feature is on. If
    feedback=`off', the feedback feature will be turned off. 
  



\section xml_page Page 

A page is just more than one panels being put together. So all its child elements
can only be panels. A page also uses attribute `numColumn'  to specify
its layout, same as panel. 

\section xml_window Specialized windows

Several special windows are designed to allow more functionality, or
to arrange panels more conveniently. These windows can be used as one
of the stacked widget in main window, or as an indpendent windows.

<ul>
 <li> \b importWindow: import window allow the user to convert other grid
    formats to volume grid format,  and to  check the quality of a
    volume grid file.  The user can modify xml files to add more grid
    types,   change default directory, or
    add usage information.
    
    To add a new grid format, one child element can be added to the
    element `gridtypes', the following attribute need/can be
    specified:

   
    - \b toXdr: if this grid type need be converted to xdr file first,
	attribute `toXdr' specifies the command.\n\n
    - \b toVog: The command that converts this grid type ( or xdr file if
	this grid type can not directly convert to .vog file) to
	.vog file\n\n
    - \b nameFilter: The name filter of this grid type.\n\n
    - \b label: The label of this grid type that appears on the interface.
        
    The arguments of the command that converts this grid type to xdr file
    is specified by child elements of the element that represents this
    grid type. Usually the format of the arguments list is
    $-variable1\ 
    value1\ -variable2\ value2 ... $, For each argument, $variablei$ is
    the tag name, the attribute `element' is used to specify the data
    type of $valuei$, the text node is used to specify the label for
    this argument.      

    The argument list  of the command that convert this grid type or .xdr
    file to .vog file is defined in the element `vogOptions'

    There is a button `Usage' on the interface, the usage information
    of a grid type can be added to the child elements of element
    `usage' as a text node.

    The default directory of grid files can be changed by modifying the attribute `dir' of
    element `directory'. 

    <li> \b physicsWindow: physicsWindow is used to define the physics state
    of the system. It interacts with other part of the system by sending signal stateChanged() and componentsChanged().
    
    have attribute element=`stateStack'. The system  state is just a
    collection of states that are specified by its child elements, so
    is its output to .vars file. In physicsWindow, the output to .vars
    file is only specified in models.  For example, Element 
    `laminar' with attribute element= `models', then all its child
    elements specifies  the models that will be written into .vars
    files if the condition is satisfied.  

    <li> \b solverWindow: solverWindow  arranges multiple
    panels by using tabbed widget.  Its child elements can only be `panel' or `page'. 
    The user need specify all panels/pages,
    the output of solverWindow to .vars file
    is a collection of  the output of all the panels/pages. 
  
    
    <li> \b initialWindow: initialWindow is for specifying compound variable \a initialConditions. 
    InitialWindow has the hard-coded part to allow the user to define regions and interact with graphics viewer. 
   




 <li>  \b boundaryWindow: boundaryWindow is for specifying
    the compound variable  \a boundary_conditions. The xml file
    only need specify all boundary types. Each boundary type is a
    panel. So all the child elements of boundaryWindow are panels.
    
    The other parts of the window is hard-coded.

    The physicsWindow will use Element 'boundary_conditions' in  xml
    file to record its current state, the user can ignore this part.
  
        

</ul>

\section xml_mainwindow Main Window

 In main window, all functionality that described in xml files are represented in a flow bar,
 which is a group of QPushButtons. Each QPushButton corresponds to a
 child element of the element `mainWindow' in xml file.  
 The child elements of mainWindow can be specialized windows, pages
 or panels. If a child element has not attribute ``element'', then its
 a hard-coded action.  The windows/panels/pages  may be displayed in central widget.
 

 Except for attribute `element',  the following attribute can be specified for the child elements of
  mainWondow:


- \b buttonTitle: the text appeared on the QPushBotton in flow bar.\n\n  
- \b title: title of the window/page/panel\n\n
- \b label: if the window/page/panel has a label on it that tell the user where to start,
    this attribute will need to be specified.\n\n
- \b status: show how much the user has done with the page. If a child
    element is required, then initialize the status to `new', otherwise
    , set the initial value to `done'. \n\n
- \b inDock: Whether graphics viewer need to be displayed. 
*/




