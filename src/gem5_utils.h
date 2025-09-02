void print_gem5_simulated_time(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return;
    }

    char line[512];

    while (fgets(line, sizeof(line), file)) {
        // Skip lines that don't start with "sim_seconds"
        if (strncmp(line, "sim_seconds", strlen("sim_seconds")) != 0)
            continue;

        // Tokenize the line using space, tab, and newline as delimiters
        char *token = strtok(line, " \t\n"); // first token
        token = strtok(NULL, " \t\n");       // second token

        if (token) {
            printf("Simulated Time : %s seconds\n", token);
        }

        break; // Only process the first matching line
    }

    fclose(file);
}

// YU modified
void print_gem5_simulated_time_part(const char *filename, int k) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return;
    }

    char line[512];
    float sec = 0;
    int cnt = 0;

    while (fgets(line, sizeof(line), file)) {
        // Skip lines that don't start with "sim_seconds"
        if (strncmp(line, "sim_seconds", strlen("sim_seconds")) != 0)
            continue;

        // Tokenize the line using space, tab, and newline as delimiters
        char *token = strtok(line, " \t\n"); // first token
        token = strtok(NULL, " \t\n");       // second token

        if (token) {
            //printf("Simulated Time : %s seconds\n", token);
            sec += atof(token);
            cnt++;
	    //printf("%d %f\n", cnt);
        }
        
        if (cnt == k+2) {
            printf("Simulated Time : %f seconds\n", sec);
             break; // Only process the first matching line
        }
           
    }

    fclose(file);
    return;
}