#include <stdio.h>
#include "node.h"
#include <string.h>
int Peek(node* top, char* title, char* composer, int* yearComposed) {
    if (top == NULL) {
        // Empty stack
        return 0;
    }

    strcpy(title, top->title);
    strcpy(composer, top->composer);
    *yearComposed = top->yearComposed;
    return 1;  // Music piece is available
}
