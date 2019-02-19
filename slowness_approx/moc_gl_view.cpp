/****************************************************************************
** Meta object code from reading C++ file 'gl_view.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "gl_view.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'gl_view.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_GlView_t {
    QByteArrayData data[16];
    char stringdata0[162];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_GlView_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_GlView_t qt_meta_stringdata_GlView = {
    {
QT_MOC_LITERAL(0, 0, 6), // "GlView"
QT_MOC_LITERAL(1, 7, 10), // "changeXang"
QT_MOC_LITERAL(2, 18, 0), // ""
QT_MOC_LITERAL(3, 19, 10), // "changeYang"
QT_MOC_LITERAL(4, 30, 10), // "changeZang"
QT_MOC_LITERAL(5, 41, 9), // "setShiftx"
QT_MOC_LITERAL(6, 51, 9), // "setShifty"
QT_MOC_LITERAL(7, 61, 11), // "setVisible0"
QT_MOC_LITERAL(8, 73, 11), // "setVisible1"
QT_MOC_LITERAL(9, 85, 11), // "setVisible2"
QT_MOC_LITERAL(10, 97, 11), // "changeTheta"
QT_MOC_LITERAL(11, 109, 9), // "changePhi"
QT_MOC_LITERAL(12, 119, 6), // "doTest"
QT_MOC_LITERAL(13, 126, 9), // "showSlice"
QT_MOC_LITERAL(14, 136, 10), // "polarSlice"
QT_MOC_LITERAL(15, 147, 14) // "setSurfaceWire"

    },
    "GlView\0changeXang\0\0changeYang\0changeZang\0"
    "setShiftx\0setShifty\0setVisible0\0"
    "setVisible1\0setVisible2\0changeTheta\0"
    "changePhi\0doTest\0showSlice\0polarSlice\0"
    "setSurfaceWire"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_GlView[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      14,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,   84,    2, 0x0a /* Public */,
       3,    1,   87,    2, 0x0a /* Public */,
       4,    1,   90,    2, 0x0a /* Public */,
       5,    1,   93,    2, 0x0a /* Public */,
       6,    1,   96,    2, 0x0a /* Public */,
       7,    1,   99,    2, 0x0a /* Public */,
       8,    1,  102,    2, 0x0a /* Public */,
       9,    1,  105,    2, 0x0a /* Public */,
      10,    1,  108,    2, 0x0a /* Public */,
      11,    1,  111,    2, 0x0a /* Public */,
      12,    0,  114,    2, 0x0a /* Public */,
      13,    0,  115,    2, 0x0a /* Public */,
      14,    0,  116,    2, 0x0a /* Public */,
      15,    1,  117,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,

       0        // eod
};

void GlView::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        GlView *_t = static_cast<GlView *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->changeXang((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->changeYang((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->changeZang((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: _t->setShiftx((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: _t->setShifty((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: _t->setVisible0((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 6: _t->setVisible1((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 7: _t->setVisible2((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 8: _t->changeTheta((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: _t->changePhi((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: _t->doTest(); break;
        case 11: _t->showSlice(); break;
        case 12: _t->polarSlice(); break;
        case 13: _t->setSurfaceWire((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObject GlView::staticMetaObject = {
    { &QGLWidget::staticMetaObject, qt_meta_stringdata_GlView.data,
      qt_meta_data_GlView,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *GlView::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *GlView::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_GlView.stringdata0))
        return static_cast<void*>(this);
    return QGLWidget::qt_metacast(_clname);
}

int GlView::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 14)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 14;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 14)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 14;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
