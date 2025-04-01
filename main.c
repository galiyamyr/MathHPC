#include "linked_list.h"

int main() {
    node* head = NULL;
    int num_nodes = 3;  

    // Generate and print the linked list
    GenerateList(&head, num_nodes);
    PrintList(head);

    char target = 'B';
    int position = SearchList(head, target);

    if (position != -1) {
        printf("The letter '%c' was found at position %d.\n", target, position);
    } else {
        printf("The letter '%c' was not found in the list.\n", target);
    }

    return 0;
}
