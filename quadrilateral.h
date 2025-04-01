#ifndef QUADRILATERAL_H
#define QUADRILATERAL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct point {
    double x;
    double y;
} point;

typedef struct quadrilateral {
    point node1;
    point node2;
    point node3;
    point node4;
    void (*quadrilateral_area)(struct quadrilateral* quad, double *area);
    void (*quadrilateral_perimeter)(struct quadrilateral* quad, double *perimeter);  // Add perimeter function pointer
    void (*quadrilateral_angles)(struct quadrilateral* quad, double *angle1, double *angle2, double *angle3, double *angle4);  // Add angles function pointer
} quadrilateral;

#endif // QUADRILATERAL_H
