#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "node.h"

// Function to query the option from the user
int QueryOption() {
    int option;
    printf("ENTER CHOICE: ");
    int result = scanf("%d", &option);
    if (result != 1) {
        printf("Invalid input. Try again.\n");
        // Clear input buffer
        while(getchar() != '\n'); 
    }
    return option;
}

// Function to execute the selected option
void ExecuteOption(const int option, node** top) {
    int yearComposed, stack_size;
    char title[100], composer[100];

    switch (option) {
        case 1: // Push new classical music piece to stack
            // Ask for music piece details
            printf("Enter music title: ");
            getchar();  // To consume leftover newline character
            fgets(title, 100, stdin);
            title[strcspn(title, "\n")] = 0; // Remove newline character at the end
            printf("Enter composer name: ");
            fgets(composer, 100, stdin);
            composer[strcspn(composer, "\n")] = 0; // Remove newline character
            printf("Enter year composed: ");
            scanf("%d", &yearComposed);
            Push(title, composer, yearComposed, top);  // Push new music piece
            break;
        
        case 2: // Pop top classical music piece off stack
            if (*top != NULL) {
                Pop(top, title, composer, &yearComposed);  // Pop the top music piece
                printf("Popped music piece: Title='%s', Composer='%s', Year Composed=%d\n", title, composer, yearComposed);
            } else {
                printf("Stack is empty.\n");
            }
            break;

        case 3: // Peek at top classical music piece
            if (*top != NULL) {
                if (Peek(*top, title, composer, &yearComposed)) {  // Peek at the top piece
                    printf("Top music piece: Title='%s', Composer='%s', Year Composed=%d\n", title, composer, yearComposed);
                }
            } else {
                printf("Stack is empty.\n");
            }
            break;

        case 4: // Display entire stack
            DisplayStack(*top);  // Display all music pieces in stack
            break;

        case 5: // Print stack size
            GetStackSize(*top, &stack_size);  // Get stack size
            printf("Stack size is %d\n", stack_size);
            break;

        case 6: // Exit (do nothing)
            break;

        default:
            printf("Error: Incorrect option. Try again.\n");
            break;
    }
}
