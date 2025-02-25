#ifndef LINKED_LIST_H
#define LINKED_LIST_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Node structure
typedef struct node {
    int value;
    struct node* next;
} node;

// Function prototypes
void GenerateList(node** head, int num_nodes);
void PrintList(node* head);

#endif // LINKED_LIST_H
