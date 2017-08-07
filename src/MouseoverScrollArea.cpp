#include "MouseoverScrollArea.hpp"

void MouseoverScrollArea::mousePressEvent(QMouseEvent *event)
{
    emit scrollAreaMouseClick((float) event->x(), (float) event->y());
}

void MouseoverScrollArea::setScale(double scale)
{
    this->scale = scale;
}
