/****************************************************************************
** Meta object code from reading C++ file 'pages.h'
**
** Created: Mon Sep 14 09:59:16 2009
**      by: The Qt Meta Object Compiler version 59 (Qt 4.4.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "pages.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'pages.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.4.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_GeneralWindow[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
      15,   14,   14,   14, 0x05,
      37,   14,   14,   14, 0x05,
      58,   14,   14,   14, 0x05,
      73,   14,   14,   14, 0x05,

 // slots: signature, parameters, type, tag, flags
      93,   14,   14,   14, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_GeneralWindow[] = {
    "GeneralWindow\0\0updateStatus(QString)\0"
    "updateStatusTip(int)\0stateChanged()\0"
    "componentsChanged()\0changeState()\0"
};

const QMetaObject GeneralWindow::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_GeneralWindow,
      qt_meta_data_GeneralWindow, 0 }
};

const QMetaObject *GeneralWindow::metaObject() const
{
    return &staticMetaObject;
}

void *GeneralWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_GeneralWindow))
        return static_cast<void*>(const_cast< GeneralWindow*>(this));
    return QWidget::qt_metacast(_clname);
}

int GeneralWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: updateStatus((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 1: updateStatusTip((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: stateChanged(); break;
        case 3: componentsChanged(); break;
        case 4: changeState(); break;
        }
        _id -= 5;
    }
    return _id;
}

// SIGNAL 0
void GeneralWindow::updateStatus(const QString & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void GeneralWindow::updateStatusTip(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void GeneralWindow::stateChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 2, 0);
}

// SIGNAL 3
void GeneralWindow::componentsChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 3, 0);
}
static const uint qt_meta_data_FloatEdit[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
      11,   10,   10,   10, 0x05,

 // slots: signature, parameters, type, tag, flags
      32,   10,   10,   10, 0x0a,
      46,   10,   10,   10, 0x0a,
      63,   10,   10,   10, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_FloatEdit[] = {
    "FloatEdit\0\0valueChanged(double)\0"
    "setValue(int)\0setValue(double)\0"
    "changeValue(QString)\0"
};

const QMetaObject FloatEdit::staticMetaObject = {
    { &QLineEdit::staticMetaObject, qt_meta_stringdata_FloatEdit,
      qt_meta_data_FloatEdit, 0 }
};

const QMetaObject *FloatEdit::metaObject() const
{
    return &staticMetaObject;
}

void *FloatEdit::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_FloatEdit))
        return static_cast<void*>(const_cast< FloatEdit*>(this));
    return QLineEdit::qt_metacast(_clname);
}

int FloatEdit::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QLineEdit::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: valueChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 1: setValue((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: setValue((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 3: changeValue((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        }
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void FloatEdit::valueChanged(double _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_IntEdit[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
       9,    8,    8,    8, 0x05,

 // slots: signature, parameters, type, tag, flags
      27,    8,    8,    8, 0x0a,
      41,    8,    8,    8, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_IntEdit[] = {
    "IntEdit\0\0valueChanged(int)\0setValue(int)\0"
    "changeValue(QString)\0"
};

const QMetaObject IntEdit::staticMetaObject = {
    { &QLineEdit::staticMetaObject, qt_meta_stringdata_IntEdit,
      qt_meta_data_IntEdit, 0 }
};

const QMetaObject *IntEdit::metaObject() const
{
    return &staticMetaObject;
}

void *IntEdit::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_IntEdit))
        return static_cast<void*>(const_cast< IntEdit*>(this));
    return QLineEdit::qt_metacast(_clname);
}

int IntEdit::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QLineEdit::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: setValue((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: changeValue((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        }
        _id -= 3;
    }
    return _id;
}

// SIGNAL 0
void IntEdit::valueChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_GeneralGroup[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
      14,   13,   13,   13, 0x05,
      29,   13,   13,   13, 0x05,
      49,   13,   13,   13, 0x05,

 // slots: signature, parameters, type, tag, flags
      70,   13,   13,   13, 0x0a,
      84,   13,   13,   13, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_GeneralGroup[] = {
    "GeneralGroup\0\0stateChanged()\0"
    "componentsChanged()\0textChanged(QString)\0"
    "changeState()\0updateChecked()\0"
};

const QMetaObject GeneralGroup::staticMetaObject = {
    { &QGroupBox::staticMetaObject, qt_meta_stringdata_GeneralGroup,
      qt_meta_data_GeneralGroup, 0 }
};

const QMetaObject *GeneralGroup::metaObject() const
{
    return &staticMetaObject;
}

void *GeneralGroup::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_GeneralGroup))
        return static_cast<void*>(const_cast< GeneralGroup*>(this));
    return QGroupBox::qt_metacast(_clname);
}

int GeneralGroup::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGroupBox::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: stateChanged(); break;
        case 1: componentsChanged(); break;
        case 2: textChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 3: changeState(); break;
        case 4: updateChecked(); break;
        }
        _id -= 5;
    }
    return _id;
}

// SIGNAL 0
void GeneralGroup::stateChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}

// SIGNAL 1
void GeneralGroup::componentsChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}

// SIGNAL 2
void GeneralGroup::textChanged(const QString & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}
static const uint qt_meta_data_OpGroup[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
      11,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
       9,    8,    8,    8, 0x0a,
      28,    8,    8,    8, 0x0a,
      42,    8,    8,    8, 0x08,
      66,    8,    8,    8, 0x08,
      90,    8,    8,    8, 0x08,
     114,    8,    8,    8, 0x08,
     141,    8,    8,    8, 0x08,
     165,    8,    8,    8, 0x08,
     183,    8,    8,    8, 0x08,
     206,    8,    8,    8, 0x08,
     237,  227,    8,    8, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_OpGroup[] = {
    "OpGroup\0\0updateComponents()\0changeState()\0"
    "updateCurrentX(QString)\0updateCurrentY(QString)\0"
    "updateCurrentZ(QString)\0"
    "updateCurrentUnit(QString)\0"
    "updateUnitList(QString)\0updateLabels(int)\0"
    "updateCurrent(QString)\0updateSelection(int)\0"
    "satisfied\0setDefault(bool)\0"
};

const QMetaObject OpGroup::staticMetaObject = {
    { &GeneralGroup::staticMetaObject, qt_meta_stringdata_OpGroup,
      qt_meta_data_OpGroup, 0 }
};

const QMetaObject *OpGroup::metaObject() const
{
    return &staticMetaObject;
}

void *OpGroup::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_OpGroup))
        return static_cast<void*>(const_cast< OpGroup*>(this));
    return GeneralGroup::qt_metacast(_clname);
}

int OpGroup::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = GeneralGroup::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: updateComponents(); break;
        case 1: changeState(); break;
        case 2: updateCurrentX((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 3: updateCurrentY((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 4: updateCurrentZ((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 5: updateCurrentUnit((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 6: updateUnitList((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 7: updateLabels((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 8: updateCurrent((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 9: updateSelection((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: setDefault((*reinterpret_cast< bool(*)>(_a[1]))); break;
        }
        _id -= 11;
    }
    return _id;
}
static const uint qt_meta_data_AllVBGroup[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      14,   12,   11,   11, 0x0a,
      35,   11,   11,   11, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_AllVBGroup[] = {
    "AllVBGroup\0\0i\0childrenClicked(int)\0"
    "updateCurrentText()\0"
};

const QMetaObject AllVBGroup::staticMetaObject = {
    { &GeneralGroup::staticMetaObject, qt_meta_stringdata_AllVBGroup,
      qt_meta_data_AllVBGroup, 0 }
};

const QMetaObject *AllVBGroup::metaObject() const
{
    return &staticMetaObject;
}

void *AllVBGroup::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_AllVBGroup))
        return static_cast<void*>(const_cast< AllVBGroup*>(this));
    return GeneralGroup::qt_metacast(_clname);
}

int AllVBGroup::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = GeneralGroup::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: childrenClicked((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: updateCurrentText(); break;
        }
        _id -= 2;
    }
    return _id;
}
static const uint qt_meta_data_StackGroup[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x0a,
      32,   11,   11,   11, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_StackGroup[] = {
    "StackGroup\0\0updateCurrentText()\0"
    "changePage(int)\0"
};

const QMetaObject StackGroup::staticMetaObject = {
    { &GeneralGroup::staticMetaObject, qt_meta_stringdata_StackGroup,
      qt_meta_data_StackGroup, 0 }
};

const QMetaObject *StackGroup::metaObject() const
{
    return &staticMetaObject;
}

void *StackGroup::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_StackGroup))
        return static_cast<void*>(const_cast< StackGroup*>(this));
    return GeneralGroup::qt_metacast(_clname);
}

int StackGroup::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = GeneralGroup::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: updateCurrentText(); break;
        case 1: changePage((*reinterpret_cast< int(*)>(_a[1]))); break;
        }
        _id -= 2;
    }
    return _id;
}
static const uint qt_meta_data_OptionPage[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x05,
      33,   11,   11,   11, 0x05,

 // slots: signature, parameters, type, tag, flags
      55,   11,   11,   11, 0x0a,
      82,   11,   11,   11, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_OptionPage[] = {
    "OptionPage\0\0updateStatusTip(int)\0"
    "updateStatus(QString)\0advancedButtonClicked(int)\0"
    "updateCurrentText()\0"
};

const QMetaObject OptionPage::staticMetaObject = {
    { &GeneralGroup::staticMetaObject, qt_meta_stringdata_OptionPage,
      qt_meta_data_OptionPage, 0 }
};

const QMetaObject *OptionPage::metaObject() const
{
    return &staticMetaObject;
}

void *OptionPage::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_OptionPage))
        return static_cast<void*>(const_cast< OptionPage*>(this));
    return GeneralGroup::qt_metacast(_clname);
}

int OptionPage::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = GeneralGroup::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: updateStatusTip((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: updateStatus((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 2: advancedButtonClicked((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: updateCurrentText(); break;
        }
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void OptionPage::updateStatusTip(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void OptionPage::updateStatus(QString _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}
static const uint qt_meta_data_AllVBWindow[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      18,   12,   13,   12, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_AllVBWindow[] = {
    "AllVBWindow\0\0bool\0save()\0"
};

const QMetaObject AllVBWindow::staticMetaObject = {
    { &GeneralWindow::staticMetaObject, qt_meta_stringdata_AllVBWindow,
      qt_meta_data_AllVBWindow, 0 }
};

const QMetaObject *AllVBWindow::metaObject() const
{
    return &staticMetaObject;
}

void *AllVBWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_AllVBWindow))
        return static_cast<void*>(const_cast< AllVBWindow*>(this));
    return GeneralWindow::qt_metacast(_clname);
}

int AllVBWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = GeneralWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: { bool _r = save();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        }
        _id -= 1;
    }
    return _id;
}
static const uint qt_meta_data_ChoiceGroup[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      13,   12,   12,   12, 0x0a,
      33,   12,   12,   12, 0x0a,
      45,   12,   12,   12, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_ChoiceGroup[] = {
    "ChoiceGroup\0\0editButtonPressed()\0"
    "update(int)\0updateCurrentText()\0"
};

const QMetaObject ChoiceGroup::staticMetaObject = {
    { &GeneralGroup::staticMetaObject, qt_meta_stringdata_ChoiceGroup,
      qt_meta_data_ChoiceGroup, 0 }
};

const QMetaObject *ChoiceGroup::metaObject() const
{
    return &staticMetaObject;
}

void *ChoiceGroup::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ChoiceGroup))
        return static_cast<void*>(const_cast< ChoiceGroup*>(this));
    return GeneralGroup::qt_metacast(_clname);
}

int ChoiceGroup::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = GeneralGroup::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: editButtonPressed(); break;
        case 1: update((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: updateCurrentText(); break;
        }
        _id -= 3;
    }
    return _id;
}
static const uint qt_meta_data_StackGroup2[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
      10,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
      13,   12,   12,   12, 0x05,

 // slots: signature, parameters, type, tag, flags
      35,   12,   12,   12, 0x0a,
      65,   63,   12,   12, 0x0a,
      81,   12,   12,   12, 0x0a,
     110,   12,  102,   12, 0x0a,
     120,   12,   12,   12, 0x0a,
     148,   12,   12,   12, 0x0a,
     154,   12,   12,   12, 0x0a,
     174,   12,   12,   12, 0x0a,
     193,   12,   12,   12, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_StackGroup2[] = {
    "StackGroup2\0\0stateChanged(QString)\0"
    "parentStateChanged(QString)\0i\0"
    "changePage(int)\0updateState(QString)\0"
    "QString\0myState()\0setParentState(QStringList)\0"
    "add()\0updateCurrentText()\0updateComponents()\0"
    "changeState()\0"
};

const QMetaObject StackGroup2::staticMetaObject = {
    { &GeneralGroup::staticMetaObject, qt_meta_stringdata_StackGroup2,
      qt_meta_data_StackGroup2, 0 }
};

const QMetaObject *StackGroup2::metaObject() const
{
    return &staticMetaObject;
}

void *StackGroup2::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_StackGroup2))
        return static_cast<void*>(const_cast< StackGroup2*>(this));
    return GeneralGroup::qt_metacast(_clname);
}

int StackGroup2::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = GeneralGroup::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: stateChanged((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 1: parentStateChanged((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 2: changePage((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: updateState((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 4: { QString _r = myState();
            if (_a[0]) *reinterpret_cast< QString*>(_a[0]) = _r; }  break;
        case 5: setParentState((*reinterpret_cast< const QStringList(*)>(_a[1]))); break;
        case 6: add(); break;
        case 7: updateCurrentText(); break;
        case 8: updateComponents(); break;
        case 9: changeState(); break;
        }
        _id -= 10;
    }
    return _id;
}

// SIGNAL 0
void StackGroup2::stateChanged(QString _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_Page[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
       6,    5,    5,    5, 0x05,

 // slots: signature, parameters, type, tag, flags
      27,    5,    5,    5, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_Page[] = {
    "Page\0\0textChanged(QString)\0"
    "updateCurrentText()\0"
};

const QMetaObject Page::staticMetaObject = {
    { &GeneralWindow::staticMetaObject, qt_meta_stringdata_Page,
      qt_meta_data_Page, 0 }
};

const QMetaObject *Page::metaObject() const
{
    return &staticMetaObject;
}

void *Page::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Page))
        return static_cast<void*>(const_cast< Page*>(this));
    return GeneralWindow::qt_metacast(_clname);
}

int Page::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = GeneralWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: textChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 1: updateCurrentText(); break;
        }
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void Page::textChanged(const QString & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
