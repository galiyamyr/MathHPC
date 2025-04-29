#include "quadrilateral.h"

// Function to calculate the angle between two vectors
double angle_between_vectors(point p1, point p2, point p3) {
    double ux = p2.x - p1.x;
    double uy = p2.y - p1.y;
    double vx = p3.x - p2.x;
    double vy = p3.y - p2.y;
    
    double dot_product = ux * vx + uy * vy;
    double magnitude_u = sqrt(ux * ux + uy * uy);
    double magnitude_v = sqrt(vx * vx + vy * vy);
    
    return acos(dot_product / (magnitude_u * magnitude_v)) * (180.0 / M_PI); // Angle in degrees
}

void quadrilateral_angles(quadrilateral* quad, double *angle1, double *angle2, double *angle3, double *angle4) {
    *angle1 = angle_between_vectors(quad->node1, quad->node2, quad->node3);
    *angle2 = angle_between_vectors(quad->node2, quad->node3, quad->node4);
    *angle3 = angle_between_vectors(quad->node3, quad->node4, quad->node1);
    *angle4 = angle_between_vectors(quad->node4, quad->node1, quad->node2);
}
