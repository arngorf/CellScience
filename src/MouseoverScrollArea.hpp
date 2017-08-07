#ifndef MOUSEOVERSCROLLAREA_HPP
#define MOUSEOVERSCROLLAREA_HPP

#include <QImage>
#include <QScrollArea>
#include <iostream>
#include <QEvent>
#include <QMouseEvent>

class MouseoverScrollArea : public QScrollArea
{
    Q_OBJECT

public:
    MouseoverScrollArea() : QScrollArea()
    {
        scale = 1.0;
    };

    void mousePressEvent(QMouseEvent *event);

    void setScale(double scale);

signals:

    void scrollAreaMouseClick(float x, float y);

private:
    double scale;
};

#endif // MOUSEOVERSCROLLAREA_HPP