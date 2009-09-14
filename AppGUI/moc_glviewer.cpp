/****************************************************************************
** Meta object code from reading C++ file 'glviewer.h'
**
** Created: Mon Sep 14 09:59:16 2009
**      by: The Qt Meta Object Compiler version 59 (Qt 4.4.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "glviewer.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'glviewer.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.4.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_GLViewer[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
      17,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
      10,    9,    9,    9, 0x05,

 // slots: signature, parameters, type, tag, flags
      27,    9,    9,    9, 0x0a,
      44,    9,    9,    9, 0x0a,
      57,    9,    9,    9, 0x0a,
      73,    9,    9,    9, 0x0a,
      89,    9,    9,    9, 0x0a,
     105,    9,    9,    9, 0x0a,
     128,  121,    9,    9, 0x0a,
     148,    9,    9,    9, 0x0a,
     154,    9,    9,    9, 0x0a,
     166,  162,    9,    9, 0x0a,
     194,  192,    9,    9, 0x0a,
     218,    9,    9,    9, 0x0a,
     239,    9,    9,    9, 0x0a,
     258,  254,    9,    9, 0x0a,
     285,    9,    9,    9, 0x0a,
     293,    9,    9,    9, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_GLViewer[] = {
    "GLViewer\0\0pickCurrent(int)\0toggleContours()\0"
    "toggleGrid()\0toggleShading()\0"
    "setShadeType1()\0setShadeType2()\0"
    "setShadeType3()\0number\0changeContours(int)\0"
    "cut()\0uncut()\0i,c\0setCurrentObj(int,QColor)\0"
    ",\0setVisibility(int,bool)\0"
    "showBoundaries(bool)\0clearCurrent()\0"
    "Nfo\0previewCut(cutplane_info&)\0reset()\0"
    "fit()\0"
};

const QMetaObject GLViewer::staticMetaObject = {
    { &QGLWidget::staticMetaObject, qt_meta_stringdata_GLViewer,
      qt_meta_data_GLViewer, 0 }
};

const QMetaObject *GLViewer::metaObject() const
{
    return &staticMetaObject;
}

void *GLViewer::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_GLViewer))
        return static_cast<void*>(const_cast< GLViewer*>(this));
    return QGLWidget::qt_metacast(_clname);
}

int GLViewer::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: pickCurrent((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: toggleContours(); break;
        case 2: toggleGrid(); break;
        case 3: toggleShading(); break;
        case 4: setShadeType1(); break;
        case 5: setShadeType2(); break;
        case 6: setShadeType3(); break;
        case 7: changeContours((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 8: cut(); break;
        case 9: uncut(); break;
        case 10: setCurrentObj((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< QColor(*)>(_a[2]))); break;
        case 11: setVisibility((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< bool(*)>(_a[2]))); break;
        case 12: showBoundaries((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 13: clearCurrent(); break;
        case 14: previewCut((*reinterpret_cast< cutplane_info(*)>(_a[1]))); break;
        case 15: reset(); break;
        case 16: fit(); break;
        }
        _id -= 17;
    }
    return _id;
}

// SIGNAL 0
void GLViewer::pickCurrent(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
