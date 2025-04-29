#include "linked_list.h"
#include <stdio.h>

// Function to search for a letter in the linked list
int SearchList(node* head, char target) {
    node* current = head;
    while (current != NULL) {
        if (current->letter == target) {  // Compare with 'letter' instead of 'value'
            return current->position;  // Return the position if found
        }
        current = current->next;
    }
    return -1;  // Return -1 if the letter is not found
}
