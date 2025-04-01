#ifndef LINKED_LIST_H
#define LINKED_LIST_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Node structure
typedef struct node {
    char letter;   // This will store the letter (instead of `value` as in the original code)
    int position;  // Position of the letter in the list
    struct node* next; // Pointer to the next node
} node;

// Function prototypes
void GenerateList(node** head, int num_nodes);  // Generates a linked list with random letters
void PrintList(node* head);  // Prints the linked list
int SearchList(node* head, char target);  // Searches for a specific letter in the linked list

#endif // LINKED_LIST_H
