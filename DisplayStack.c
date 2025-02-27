#include <stdio.h>
#include "node.h"

void DisplayStack(node* top) {
    if (top == NULL) {
        printf("Stack is empty.\n");
        return;
    }

    printf("Music Pieces in the stack:\n");
    while (top != NULL) {
        printf("Title: '%s', Composer: '%s', Year Composed: %d, Position: %d\n",
               top->title, top->composer, top->yearComposed, top->position);
        top = top->next;
    }
}
