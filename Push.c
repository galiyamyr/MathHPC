#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "node.h"

void Push(const char* title, const char* composer, const int yearComposed, node** top) {
    node* new_node = (node*)malloc(sizeof(node));
    strcpy(new_node->title, title);
    strcpy(new_node->composer, composer);
    new_node->yearComposed = yearComposed;
    new_node->position = (*top == NULL) ? 1 : (*top)->position + 1;
    new_node->next = *top;
    *top = new_node;
}
