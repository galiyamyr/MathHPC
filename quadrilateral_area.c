#include "quadrilateral.h"

double triangle_area(point p1, point p2, point p3) {
    return fabs((p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y)) / 2.0);
}

void quadrilateral_area(quadrilateral* quad, double *area) {
    double area1 = triangle_area(quad->node1, quad->node2, quad->node3);
    double area2 = triangle_area(quad->node1, quad->node3, quad->node4);
    *area = area1 + area2;
}
