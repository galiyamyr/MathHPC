#include "quadrilateral.h"

// Function to calculate distance between two points
double distance(point p1, point p2) {
    return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

void quadrilateral_perimeter(quadrilateral* quad, double *perimeter) {
    double side1 = distance(quad->node1, quad->node2);
    double side2 = distance(quad->node2, quad->node3);
    double side3 = distance(quad->node3, quad->node4);
    double side4 = distance(quad->node4, quad->node1);
    
    *perimeter = side1 + side2 + side3 + side4;
}
