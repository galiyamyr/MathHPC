#include <stdio.h>
#include <stdlib.h>
#include "node.h"  

int main() {
    node* top = NULL;  // Initialize the stack as empty
    int option;

    // Display options first
    DisplayOptions();

    // Loop until user chooses to exit
    while (1) {
        // Query user for an option
        option = QueryOption();

        // If the user selects 6, exit the loop
        if (option == 6) {
            printf("Exiting program...\n");
            break;
        }

        // Execute the chosen option
        ExecuteOption(option, &top);

        // Optionally, re-display the menu after each action
        DisplayOptions();
    }

    // Clean up stack before exiting
    DeleteStack(&top);

    return 0;
}
