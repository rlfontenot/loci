/****************************************************************************
** Meta object code from reading C++ file 'cutdialog.h'
**
** Created: Mon Sep 14 09:59:17 2009
**      by: The Qt Meta Object Compiler version 59 (Qt 4.4.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "cutdialog.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'cutdialog.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.4.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_CutDialog[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
      11,   10,   10,   10, 0x05,
      42,   10,   10,   10, 0x05,

 // slots: signature, parameters, type, tag, flags
      55,   10,   10,   10, 0x08,
      65,   10,   10,   10, 0x08,
      84,   10,   10,   10, 0x08,
      92,   10,   10,   10, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_CutDialog[] = {
    "CutDialog\0\0cutInfoChanged(cutplane_info&)\0"
    "cutPressed()\0setInfo()\0planeSelected(int)\0"
    "reset()\0cut()\0"
};

const QMetaObject CutDialog::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_CutDialog,
      qt_meta_data_CutDialog, 0 }
};

const QMetaObject *CutDialog::metaObject() const
{
    return &staticMetaObject;
}

void *CutDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_CutDialog))
        return static_cast<void*>(const_cast< CutDialog*>(this));
    return QWidget::qt_metacast(_clname);
}

int CutDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: cutInfoChanged((*reinterpret_cast< cutplane_info(*)>(_a[1]))); break;
        case 1: cutPressed(); break;
        case 2: setInfo(); break;
        case 3: planeSelected((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: reset(); break;
        case 5: cut(); break;
        }
        _id -= 6;
    }
    return _id;
}

// SIGNAL 0
void CutDialog::cutInfoChanged(cutplane_info & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void CutDialog::cutPressed()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}
QT_END_MOC_NAMESPACE
