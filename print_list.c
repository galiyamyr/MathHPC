#include <stdio.h>
#include <stdlib.h>
#include "linked_list.h"

void PrintList(node* head) {
    node* temp = head;
    int pos = 1;

    // Print header
    printf("\n----------------------------------------------------\n");
    printf("| Pos | Val |    Address    |     Next      |\n");
    printf("----------------------------------------------------\n");

    // Iterate through the linked list
    while (temp) {
        printf("| %-3d | %-3c | %p | %p |\n", pos, temp->letter, (void*)temp, (void*)temp->next);
        temp = temp->next;
        pos++;
    }

    printf("----------------------------------------------------\n");
}
