#include <stdio.h>
#include <stdlib.h>
#include "node.h"
#include <string.h>
void Pop(node** top, char* title, char* composer, int* yearComposed) {
    if (*top == NULL) {
        // Indicating an empty stack
        return;
    }

    node* temp = *top;
    strcpy(title, temp->title);
    strcpy(composer, temp->composer);
    *yearComposed = temp->yearComposed;

    *top = temp->next;
    free(temp);
}
