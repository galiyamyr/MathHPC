#include "linked_list.h"

// Function to generate a random letter (uppercase or lowercase)
char random_letter() {
    return (rand() % 2 == 0) ? ('A' + rand() % 26) : ('a' + rand() % 26);
}

// Function to create a linked list with random letters
void GenerateList(node** head, int num_nodes) {
    if (num_nodes <= 0) return;

    srand(time(NULL));  // Seed the random number generator

    for (int i = 0; i < num_nodes; i++) {
        node* temp = (node*)malloc(sizeof(node));
        if (!temp) {
            printf("Memory allocation failed\n");
            return;
        }

        temp->letter = random_letter();  // Assign a random letter
        temp->position = i;  // Assign sequential position
        temp->next = *head;  // Insert at the beginning
        *head = temp;  // Update head pointer
    }
}
