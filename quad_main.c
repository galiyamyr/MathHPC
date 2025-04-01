#include "quadrilateral.h"

// Function prototypes
void quadrilateral_area(quadrilateral* quad, double *area);
void quadrilateral_perimeter(quadrilateral* quad, double *perimeter);
void quadrilateral_angles(quadrilateral* quad, double *angle1, double *angle2, double *angle3, double *angle4);

int main() {
    point node1, node2, node3, node4;
    
    // User input for the coordinates of four vertices
    printf("Enter coordinates of first point (x1 y1): ");
    scanf("%lf %lf", &node1.x, &node1.y);
    printf("Enter coordinates of second point (x2 y2): ");
    scanf("%lf %lf", &node2.x, &node2.y);
    printf("Enter coordinates of third point (x3 y3): ");
    scanf("%lf %lf", &node3.x, &node3.y);
    printf("Enter coordinates of fourth point (x4 y4): ");
    scanf("%lf %lf", &node4.x, &node4.y);
    
    // Create a quadrilateral struct and assign points
    quadrilateral quad = {node1, node2, node3, node4};

    // Assign function pointers
    quad.quadrilateral_area = quadrilateral_area;
    quad.quadrilateral_perimeter = quadrilateral_perimeter;
    quad.quadrilateral_angles = quadrilateral_angles;

    // Variables to hold the area, perimeter, and angles
    double total_area = 0.0;
    double perimeter = 0.0;
    double angle1, angle2, angle3, angle4;
    
    // Calculate area, perimeter, and angles using their respective function pointers
    quad.quadrilateral_area(&quad, &total_area);
    quad.quadrilateral_perimeter(&quad, &perimeter);
    quad.quadrilateral_angles(&quad, &angle1, &angle2, &angle3, &angle4);
    
    // Output the results
    printf("The area of the quadrilateral is: %.2lf\n", total_area);
    printf("The perimeter of the quadrilateral is: %.2lf\n", perimeter);
    printf("The angles of the quadrilateral are: %.2lf, %.2lf, %.2lf, %.2lf\n", angle1, angle2, angle3, angle4);
    
    return 0;
}
 