#define DEBUG

#ifndef DEBUG
    #include <gplib.h>
    #define _PPScriptMode_		// for enum mode, replace this with #define _EnumMode_
    #include "../../Include/pp_proj.h"
#else
    #include <stdio.h>
    #include <stdlib.h>
    #include <string.h>
#endif

#include "../thirdparty/cJSON/cJSON.h"

// Define matrix size (for simplicity, assuming a maximum of 4 DoFs)
#define MAX_DOF 4

// Define state-space matrices
double A[2 * MAX_DOF][2 * MAX_DOF] = {0};
double B[2 * MAX_DOF][1] = {0};
double C[1][2 * MAX_DOF] = {0};
double D[1][1] = {0};

// Function to read the entire contents of a file into a string
char *read_file(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("fopen");
        return NULL;
    }

    // Get the file size
    fseek(file, 0, SEEK_END);
    long length = ftell(file);
    fseek(file, 0, SEEK_SET);

    if (length < 0) {
        perror("ftell");
        fclose(file);
        return NULL;
    }

    // Allocate memory for the file contents
    char *buffer = (char *)malloc(length + 1);
    if (!buffer) {
        perror("malloc");
        fclose(file);
        return NULL;
    }

    // Read the file contents into the buffer
    fread(buffer, 1, length, file);
    if (ferror(file)) {
        perror("fread");
        fclose(file);
        free(buffer);
        return NULL;
    }

    buffer[length] = '\0'; // Null-terminate the string

    fclose(file);
    return buffer;
}

void buildMatrixA(cJSON *json_type){

        // DoF-wise mode: X^T=[x_1, xdot_1, x_2, xdot_2, ...]
        if (strcmp(id_source, "ground") == 0 && strcmp(id_dest, "ground") != 0) {
            // Update A matrix for ground connection
            A[2 * (atoi(id_dest) - 1)][2 * (atoi(id_dest) - 1) + 1] = 1;
            A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_dest) - 1)] = -stiffness / mass[atoi(id_dest)];
            A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_dest) - 1) + 1] = -damping / mass[atoi(id_dest)];

        } else if (strcmp(id_source, "ground") == 1 && strcmp(id_dest, "ground") != 0) {
            // Update A matrix for ground connection
            A[2 * (atoi(id_source) - 1)][2 * (atoi(id_source) - 1) + 1] = 1;
            A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_source) - 1)]     = -stiffness / mass[atoi(id_source)];
            A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_source) - 1) + 1] = -damping / mass[atoi(id_source)];

        } else if (atoi(id_source) > 0        && atoi(id_dest) > 0
                && atoi(id_source) <= MAX_DOF && atoi(id_dest) <= MAX_DOF)  {
            // Update A matrix for DoF connection
            A[2 * (atoi(id_dest) - 1)][2 * (atoi(id_dest) - 1) + 1] = 1;
            A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_dest) - 1)] -= stiffness / mass[atoi(id_dest)];
            A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_dest) - 1) + 1] -= damping / mass[atoi(id_dest)];
            A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_source) - 1)] += stiffness / mass[atoi(id_dest)];
            A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_source) - 1) + 1] += damping / mass[atoi(id_dest)];

            A[2 * (atoi(id_source) - 1)][2 * (atoi(id_source) - 1) + 1] = 1;
            A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_source) - 1)]     -= stiffness / mass[atoi(id_source)];
            A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_source) - 1) + 1] -= damping / mass[atoi(id_source)];
            A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_dest) - 1)]     += stiffness / mass[atoi(id_source)];
            A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_dest) - 1) + 1] += damping / mass[atoi(id_source)];
        }


        /////////////////////////////////////////////////////////////////////////////
        // Differential Order-wise mode: X^T=[x_1, x_2, ... , xdot_1, xdot_2, ...]
        // if (strcmp(id_source, "ground") == 0 && strcmp(id_dest, "ground") != 0) {
        //     // Update A matrix for ground connection
        //     A[(atoi(id_dest) - 1)][(atoi(id_dest) + DoF_count - 1)] = 1;
        //     A[(atoi(id_dest) + DoF_count - 1)][(atoi(id_dest) - 1)] = -stiffness / mass[atoi(id_dest)];
        //     A[(atoi(id_dest) + DoF_count - 1)][(atoi(id_dest) + DoF_count - 1)] = -damping / mass[atoi(id_dest)];

        // } else if (strcmp(id_source, "ground") == 1 && strcmp(id_dest, "ground") != 0) {
        //     // Update A matrix for ground connection
        //     A[(atoi(id_source) - 1)][(atoi(id_source) + DoF_count - 1)] = 1;
        //     A[(atoi(id_source) + DoF_count  - 1)][(atoi(id_source) - 1)]     = -stiffness / mass[atoi(id_source)];
        //     A[(atoi(id_source) + DoF_count  - 1)][(atoi(id_source) + DoF_count - 1)] = -damping / mass[atoi(id_source)];

        // } else if (atoi(id_source) > 0        && atoi(id_dest) > 0
        //         && atoi(id_source) <= MAX_DOF && atoi(id_dest) <= MAX_DOF)  {
        //     // Update A matrix for DoF connection
        //     A[(atoi(id_dest) - 1)][(atoi(id_dest) + DoF_count - 1)] = 1;
        //     A[(atoi(id_dest) + DoF_count - 1)][(atoi(id_dest) - 1)] -= stiffness / mass[atoi(id_dest)];
        //     A[(atoi(id_dest) + DoF_count - 1)][(atoi(id_dest) + DoF_count - 1)] -= damping / mass[atoi(id_dest)];
        //     A[(atoi(id_dest) + DoF_count - 1)][(atoi(id_source) - 1)] += stiffness / mass[atoi(id_dest)];
        //     A[(atoi(id_dest) + DoF_count - 1)][(atoi(id_source) + DoF_count - 1)] += damping / mass[atoi(id_dest)];

        //     A[(atoi(id_source) - 1)][(atoi(id_source) + DoF_count - 1)] = 1;
        //     A[(atoi(id_source) + DoF_count  - 1)][(atoi(id_source) - 1)]     -= stiffness / mass[atoi(id_source)];
        //     A[(atoi(id_source) + DoF_count  - 1)][(atoi(id_source) + DoF_count - 1)] -= damping / mass[atoi(id_source)];
        //     A[(atoi(id_source) + DoF_count  - 1)][(atoi(id_dest) - 1)]     += stiffness / mass[atoi(id_source)];
        //     A[(atoi(id_source) + DoF_count  - 1)][(atoi(id_dest) + DoF_count - 1)] += damping / mass[atoi(id_source)];
        // }
        /////////////////////////////////////////////////////////////////////////////
    }

}

void buildMatrix(cJSON *json_type) {
    // TODO: Add another argument to choose in which
    // format the matrix will be built

    // Extract DoFs - TODO: move to another function
    cJSON *DoFs = cJSON_GetObjectItem(json_type, "DoFs");
    cJSON *connections = cJSON_GetObjectItem(json_type, "connections");

    unsigned int DoF_count = cJSON_GetArraySize(DoFs);
    size_t connections_count = cJSON_GetArraySize(connections);

    int id;
    int idx;
    const char *id_source;
    const char *id_dest;
    double mass[DoF_count], stiffness, damping;

    if (DoFs == NULL) {
        printf("Error: 'DoFs' not found\n");
        return;
    } else if (!cJSON_IsArray(DoFs)) {
        printf("Error: 'DoFs' is not an array\n");
        return;
    }

    if (connections == NULL) {
        printf("Error: 'connections' not found\n");
        return;
    } else if (!cJSON_IsArray(connections)) {
        printf("Error: 'connections' is not an array\n");
        return;
    }


    // Initialize mass array
    for (idx = 0; idx < DoF_count; idx++) {
        cJSON *DoF = cJSON_GetArrayItem(DoFs, idx);
        id = cJSON_GetObjectItem(DoF, "id")->valueint;
        // TO DO: Validate id
        mass[id] = cJSON_GetObjectItem(DoF, "mass")->valuedouble;
    }

    for (i = 0; i < connections_count; i++) {
        // Extract connections
        cJSON *connection = cJSON_GetArrayItem(connections, i);
        id_source = cJSON_GetObjectItem(connection, "id_source")->valuestring;
        id_dest   = cJSON_GetObjectItem(connection, "id_dest")->valuestring;
        stiffness = cJSON_GetObjectItem(connection, "stiffness")->valuedouble;
        damping = cJSON_GetObjectItem(connection, "damping")->valuedouble;

        // DoF-wise mode: X^T=[x_1, xdot_1, x_2, xdot_2, ...]
        if (strcmp(id_source, "ground") == 0 && strcmp(id_dest, "ground") != 0) {
            // Update A matrix for ground connection
            A[2 * (atoi(id_dest) - 1)][2 * (atoi(id_dest) - 1) + 1] = 1;
            A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_dest) - 1)] = -stiffness / mass[atoi(id_dest)];
            A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_dest) - 1) + 1] = -damping / mass[atoi(id_dest)];

        } else if (strcmp(id_source, "ground") == 1 && strcmp(id_dest, "ground") != 0) {
            // Update A matrix for ground connection
            A[2 * (atoi(id_source) - 1)][2 * (atoi(id_source) - 1) + 1] = 1;
            A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_source) - 1)]     = -stiffness / mass[atoi(id_source)];
            A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_source) - 1) + 1] = -damping / mass[atoi(id_source)];

        } else if (atoi(id_source) > 0        && atoi(id_dest) > 0
                && atoi(id_source) <= MAX_DOF && atoi(id_dest) <= MAX_DOF)  {
            // Update A matrix for DoF connection
            A[2 * (atoi(id_dest) - 1)][2 * (atoi(id_dest) - 1) + 1] = 1;
            A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_dest) - 1)] -= stiffness / mass[atoi(id_dest)];
            A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_dest) - 1) + 1] -= damping / mass[atoi(id_dest)];
            A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_source) - 1)] += stiffness / mass[atoi(id_dest)];
            A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_source) - 1) + 1] += damping / mass[atoi(id_dest)];

            A[2 * (atoi(id_source) - 1)][2 * (atoi(id_source) - 1) + 1] = 1;
            A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_source) - 1)]     -= stiffness / mass[atoi(id_source)];
            A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_source) - 1) + 1] -= damping / mass[atoi(id_source)];
            A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_dest) - 1)]     += stiffness / mass[atoi(id_source)];
            A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_dest) - 1) + 1] += damping / mass[atoi(id_source)];
        }

}
}

// void buildMatrixA(cJSON *json_type) {



//     cJSON *connections = cJSON_GetObjectItem(json_type, "connections");
//     // if (connections != NULL) {

//     // }
//     // if (!cJSON_IsArray(connections)) {

//     // }

//     for (i = 0; i < connections_count; i++) {
//         // Extract connections
//         cJSON *connection = cJSON_GetArrayItem(connections, i);
//         id_source = cJSON_GetObjectItem(connection, "id_source")->valuestring;
//         id_dest   = cJSON_GetObjectItem(connection, "id_dest")->valuestring;
//         stiffness = cJSON_GetObjectItem(connection, "stiffness")->valuedouble;
//         damping = cJSON_GetObjectItem(connection, "damping")->valuedouble;

//         // DoF-wise mode: X^T=[x_1, xdot_1, x_2, xdot_2, ...]
//         if (strcmp(id_source, "ground") == 0 && strcmp(id_dest, "ground") != 0) {
//             // Update A matrix for ground connection
//             A[2 * (atoi(id_dest) - 1)][2 * (atoi(id_dest) - 1) + 1] = 1;
//             A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_dest) - 1)] = -stiffness / mass[atoi(id_dest)];
//             A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_dest) - 1) + 1] = -damping / mass[atoi(id_dest)];

//         } else if (strcmp(id_source, "ground") == 1 && strcmp(id_dest, "ground") != 0) {
//             // Update A matrix for ground connection
//             A[2 * (atoi(id_source) - 1)][2 * (atoi(id_source) - 1) + 1] = 1;
//             A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_source) - 1)]     = -stiffness / mass[atoi(id_source)];
//             A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_source) - 1) + 1] = -damping / mass[atoi(id_source)];

//         } else if (atoi(id_source) > 0        && atoi(id_dest) > 0
//                 && atoi(id_source) <= MAX_DOF && atoi(id_dest) <= MAX_DOF)  {
//             // Update A matrix for DoF connection
//             A[2 * (atoi(id_dest) - 1)][2 * (atoi(id_dest) - 1) + 1] = 1;
//             A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_dest) - 1)] -= stiffness / mass[atoi(id_dest)];
//             A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_dest) - 1) + 1] -= damping / mass[atoi(id_dest)];
//             A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_source) - 1)] += stiffness / mass[atoi(id_dest)];
//             A[2 * (atoi(id_dest) - 1) + 1][2 * (atoi(id_source) - 1) + 1] += damping / mass[atoi(id_dest)];

//             A[2 * (atoi(id_source) - 1)][2 * (atoi(id_source) - 1) + 1] = 1;
//             A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_source) - 1)]     -= stiffness / mass[atoi(id_source)];
//             A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_source) - 1) + 1] -= damping / mass[atoi(id_source)];
//             A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_dest) - 1)]     += stiffness / mass[atoi(id_source)];
//             A[2 * (atoi(id_source) - 1) + 1][2 * (atoi(id_dest) - 1) + 1] += damping / mass[atoi(id_source)];
//         }


//         /////////////////////////////////////////////////////////////////////////////
//         // Differential Order-wise mode: X^T=[x_1, x_2, ... , xdot_1, xdot_2, ...]
//         // if (strcmp(id_source, "ground") == 0 && strcmp(id_dest, "ground") != 0) {
//         //     // Update A matrix for ground connection
//         //     A[(atoi(id_dest) - 1)][(atoi(id_dest) + DoF_count - 1)] = 1;
//         //     A[(atoi(id_dest) + DoF_count - 1)][(atoi(id_dest) - 1)] = -stiffness / mass[atoi(id_dest)];
//         //     A[(atoi(id_dest) + DoF_count - 1)][(atoi(id_dest) + DoF_count - 1)] = -damping / mass[atoi(id_dest)];

//         // } else if (strcmp(id_source, "ground") == 1 && strcmp(id_dest, "ground") != 0) {
//         //     // Update A matrix for ground connection
//         //     A[(atoi(id_source) - 1)][(atoi(id_source) + DoF_count - 1)] = 1;
//         //     A[(atoi(id_source) + DoF_count  - 1)][(atoi(id_source) - 1)]     = -stiffness / mass[atoi(id_source)];
//         //     A[(atoi(id_source) + DoF_count  - 1)][(atoi(id_source) + DoF_count - 1)] = -damping / mass[atoi(id_source)];

//         // } else if (atoi(id_source) > 0        && atoi(id_dest) > 0
//         //         && atoi(id_source) <= MAX_DOF && atoi(id_dest) <= MAX_DOF)  {
//         //     // Update A matrix for DoF connection
//         //     A[(atoi(id_dest) - 1)][(atoi(id_dest) + DoF_count - 1)] = 1;
//         //     A[(atoi(id_dest) + DoF_count - 1)][(atoi(id_dest) - 1)] -= stiffness / mass[atoi(id_dest)];
//         //     A[(atoi(id_dest) + DoF_count - 1)][(atoi(id_dest) + DoF_count - 1)] -= damping / mass[atoi(id_dest)];
//         //     A[(atoi(id_dest) + DoF_count - 1)][(atoi(id_source) - 1)] += stiffness / mass[atoi(id_dest)];
//         //     A[(atoi(id_dest) + DoF_count - 1)][(atoi(id_source) + DoF_count - 1)] += damping / mass[atoi(id_dest)];

//         //     A[(atoi(id_source) - 1)][(atoi(id_source) + DoF_count - 1)] = 1;
//         //     A[(atoi(id_source) + DoF_count  - 1)][(atoi(id_source) - 1)]     -= stiffness / mass[atoi(id_source)];
//         //     A[(atoi(id_source) + DoF_count  - 1)][(atoi(id_source) + DoF_count - 1)] -= damping / mass[atoi(id_source)];
//         //     A[(atoi(id_source) + DoF_count  - 1)][(atoi(id_dest) - 1)]     += stiffness / mass[atoi(id_source)];
//         //     A[(atoi(id_source) + DoF_count  - 1)][(atoi(id_dest) + DoF_count - 1)] += damping / mass[atoi(id_source)];
//         // }
//         /////////////////////////////////////////////////////////////////////////////
//     }

// }

int parse_matrixA(cJSON *json_type) {
    int num_rows, num_columns;
    int i, j;

    cJSON *jsonA = cJSON_GetObjectItem(json_type, "A");
    if (jsonA == NULL) {
        printf("Error: 'A' not found\n");
        return 1;
    } else if (!cJSON_IsArray(jsonA)) {
        printf("Error: 'A' is not an array\n");
        return 1;
    }

    num_rows = cJSON_GetArraySize(jsonA);
    if (num_rows < 1 || num_rows > MAX_DOF){
        printf("Error: invalid number of rows\n");
        return 1;
    }
    num_columns = cJSON_GetArraySize(cJSON_GetArrayItem(jsonA, 0));
    if (num_columns < 1 || num_columns > MAX_DOF || num_columns != num_rows){
        printf("Error: invalid number of columns\n");
        return 1;
    }
    if (num_columns != num_rows){
        printf("Error: matrix A is not square\n");
        return 1;
    }

    for (i = 0; i < num_rows; i++) {
        cJSON *row = cJSON_GetArrayItem(jsonA, i);
        if (cJSON_GetArraySize(row) != num_columns) {
            printf("Error: row %d has the wrong number of columns\n", row->valueint);
            return 1;
        } else for (j = 0; j < num_columns; j++) {
            A[i][j]=cJSON_GetArrayItem(row, j)->valuedouble;
        }
    }
    return 0;
}

void test_print_matrix(double *matrix, int rows, int cols) {
    int i = 0;
    for (; i < rows; i++) {
        int j = 0;
        for (; j < cols; j++) {
            printf("%6.5f ", *(matrix+i*cols+j));
        }
        printf("\n");
    }
}

void handle_json_type( cJSON *json_type) {
    if (json_type == NULL) {
        printf("Error reading type\n");
        return;
    }
    if (!cJSON_IsString(json_type)) {
        printf("Error: 'json_type' is not a string\n");
        return;
    }

    if(strcmp(json_type->valuestring, "matrix") == 0) {
        parse_matrix(json_type);
    } else if(strcmp(json_type->valuestring, "2ndOrderSystem") == 0) {
        buildMatrix(json_type);
    } else if(strcmp(json_type->valuestring, "configTrajectory") == 0) {
        //PLACEHOLDER
    } else if(strcmp(json_type->valuestring, "configGather") == 0) {
        //PLACEHOLDER
    } else {
        printf("Error: invalid type\n");
    }
}


void parse_json(const char *json_string) {
    cJSON *json = cJSON_Parse(json_string);
    if (json == NULL) {
        const char *error_ptr = cJSON_GetErrorPtr();
        printf("Error parsing JSON\n");
        if (error_ptr != NULL) {
            printf("%s\n", error_ptr);
        }
        return;
    }

    cJSON *type = cJSON_GetObjectItem(json, "type");
    handle_json_type(type);

    cJSON_Delete(type);
    cJSON_Delete(json);
}

int main(void)
{
    // InitLibrary();  // Required for accessing Power PMAC library

    // // Set the environment variable
    // if (setenv("MY_PROGRAM_PATH", "/var/ftp/usrflash/Project/C Language/Background Programs/parse_json", 1) != 0) {
    //     perror("setenv");
    //     return 1;
    // }

    // char *json_string = read_file("/var/ftp/usrflash/Project/Documentation/system.json");
    char *json_string = read_file("../example/2ndOrderSystem.json");
    // char *json_string = read_file("../example/matrix.json");

    if (json_string == NULL) {
        return 1;
    }

    // Get JSON type
    // if 'MSD', then build matrices
    // TODO: if 'RLC', then build matrices
    // if 'matrices' then parse matrices

    parse_json(json_string);

    // Free the allocated memory
    free(json_string);

	// CloseLibrary();
	return 0;
}

