#include <stdlib.h>
#include "node.h"

void GetStackSize(node* top, int* stack_size) {
    if (top == NULL) {
        *stack_size = 0;
        return;
    }

    *stack_size = top->position;  // The top node holds the correct position (size)
}
