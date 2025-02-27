#ifndef __NODE_H__
#define __NODE_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // For strcpy, fgets, etc.
#include "node.h"
// Define a structure for Classical Music Piece without ID
typedef struct node {
    char title[100];       // Title of the music piece
    char composer[100];    // Composer of the music piece
    int yearComposed;      // Year the piece was composed
    int position;          // Position in the stack
    struct node* next;     // Pointer to the next music piece in the stack
} node;

// Function prototypes
void DisplayOptions();
int QueryOption();
void ExecuteOption(const int option, node** top);
void Push(const char* title, const char* composer, const int yearComposed, node** top);
void Pop(node** top, char* title, char* composer, int* yearComposed);
int Peek(node* top, char* title, char* composer, int* yearComposed);
void DisplayStack(node* top);
void GetStackSize(node* top, int* stack_size);
void DeleteStack(node** top);

#endif // __NODE_H__
