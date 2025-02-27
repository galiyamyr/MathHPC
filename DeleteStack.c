#include <stdio.h>
#include <stdlib.h>
#include "node.h"

void DeleteStack(node** top) {
    node* temp;
    while (*top != NULL) {
        temp = *top;
        *top = (*top)->next;
        free(temp);
    }

    printf("\nGoodbye!\n\n");
}
