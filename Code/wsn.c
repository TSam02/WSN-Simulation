#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#define _POSIX //only for Windows

#include <netdb.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#define MIN_PORT_AVAILABLE 2
#define MAX_ITERATION 10
#define DEFAULT_INTERVAL 5
#define DEFAULT_NUM_PORTS 5

#define MAX_NEIGHBOURS 4
#define MAX_NEARBY_NEIGH 16
#define MAX_SARRAY_SIZE 5
#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1

#define MSG_REQUEST_USAGE 0
#define MSG_SENT_USAGE 1
#define MSG_ALERT 2
#define MSG_RESOLUTION 3
#define MSG_EXIT 4
#define MSG_CART_DATA 5



/* ------------------------------------- */
/* ----- PRE-DECLARE OF FUNCTIONS ------ */
/* ------------------------------------- */
int base_station(MPI_Comm world_comm, MPI_Comm comm, int interval, int nports);
int charging_node(MPI_Comm world_comm, MPI_Comm comm, int nrows, int ncols, int interval, int nports, int *dims);

char* get_current_datetime();
int rand_update_availabilty();

void *ThreadFuncBaseRecv(void *pArg);
void *ThreadFuncUpdatePort(void *pArg);

/* ------------------------------------- */
/* ----- DECLARATION OF NEW STRUCT ----- */
/* ------------------------------------- */
typedef struct ports_record_struct {
    char *rec_timestamp;
    int total_free_port;  
} ports_record;

typedef struct node_thread_data_struct {
    int node_rank;
    int thread_rank; 
    int *curr_record_idx;
    int *port_shut_down;
} node_thread_data;

typedef struct alert_data_struct {
    char detected_timestamp[20];
    int num_msg_exchanged;
    double sending_time;
    double node_start_time;
    int num_free_port;
    int adjacent_nodes[MAX_NEIGHBOURS];
    int nearby_nodes[MAX_NEARBY_NEIGH];
    int my_node_rank;
    int adjacent_node_port_usage[MAX_NEIGHBOURS];
    char IPv4_add[16];
    char neighbour_IPv4[MAX_NEIGHBOURS][16];
} alert_data;

typedef struct base_thread_data_struct {
    double *report_recv_time;
    int *report_node_arr_world;
    int *report_node_arr_cart;
    alert_data *report_detail_arr;
    int *num_report_recv;
    int *shut_down;
} base_thread_data;

/* ------------------------------------- */
/* ---------- GLOBAL VARIABLE ---------- */
/* ------------------------------------- */
ports_record *ports_record_arr[MAX_SARRAY_SIZE];      // shared array for charging node

pthread_mutex_t shared_arr_Mutex;   // mutex for the shared array
pthread_cond_t shared_arr_Cond;     // conditional variable for shared array   

MPI_Datatype AlertType;


/* ------------------------------------- */
/* -------------- MAIN ----------------- */ 
/* ------------------------------------- */

int main(int argc, char *argv[]) {

    /* SECTION 1 - SETTING UP MPI  */
    int provided;     // Store provided thread support level
    int rank, size;   // rank and size of the mpi process
    int nrows, ncols, interval, nports;
    int dims[2];  
    MPI_Comm new_comm;

    // Initialize MPI environment
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_split(MPI_COMM_WORLD, rank == 0, 0, &new_comm);     // split the processes into 2 communicators

    /* process command line arguments*/
    if (argc == 5) {     
        nrows = atoi(argv[1]);      /* 1st argument - number of rows */
        ncols = atoi(argv[2]);      /* 2nd argument - number of columns */
        interval = atoi(argv[3]);   /* 3rd argument - interval of each iteration */
        nports = atoi(argv[4]);     /* 4th argument - number of ports in each charging node */

        dims[0] = nrows; 
        dims[1] = ncols; 
        interval = 5;
        nports = 5;

        /* Error Handling */
        if ((nrows * ncols) != size-1) {
            if (rank == 0) 
                printf("ERROR: nrows*ncols=%d * %d = %d != %d\n", nrows, ncols, nrows * ncols, size-1);
            MPI_Finalize();
            return 0;
        
        if (interval < 3) {
            if (rank == 0) 
                printf("ERROR: interval per iteration should be at least 3\n");
            MPI_Finalize();
            return 0;
        }

        if (nports < 3) {
            if (rank == 0) 
                printf("ERROR: number of ports per node should be at least 3\n");
            MPI_Finalize();
            return 0;
        }

        }
    } else {
        if (rank == 0) {
            printf("Insufficient input argumentns! Running default config.\n");
        }
        nrows = ncols = (int)sqrt(size-1);
        dims[0]=dims[1]=0;
        interval = DEFAULT_INTERVAL;
        nports = DEFAULT_NUM_PORTS;
    }

    /* SECTION 2 - DEFINING NEW DATATYPE FOR MPI */
    alert_data alert;
    MPI_Datatype type[11] = { MPI_CHAR, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_CHAR, MPI_CHAR };
    int blocklen[11] = {20, 1, 1, 1, 1, 4, 16, 1, 4, 16, 64};
    MPI_Aint disp[11];

    MPI_Get_address(&alert.detected_timestamp, &disp[0]);
	MPI_Get_address(&alert.num_msg_exchanged, &disp[1]);
    MPI_Get_address(&alert.sending_time, &disp[2]);
    MPI_Get_address(&alert.node_start_time, &disp[3]);
    MPI_Get_address(&alert.num_free_port, &disp[4]);
    MPI_Get_address(&alert.adjacent_nodes, &disp[5]);
    MPI_Get_address(&alert.nearby_nodes, &disp[6]);
    MPI_Get_address(&alert.my_node_rank, &disp[7]);
    MPI_Get_address(&alert.adjacent_node_port_usage, &disp[8]);
    MPI_Get_address(&alert.IPv4_add, &disp[9]);
    MPI_Get_address(&alert.neighbour_IPv4, &disp[10]);


    //Make relative
    for (int i=10; i >= 0; i--) {
        disp[i] -= disp[0];
    }
	disp[0]=0;
    
    // Create MPI struct
	MPI_Type_create_struct(11, blocklen, disp, type, &AlertType);
	MPI_Type_commit(&AlertType);


    /* SECTION 3 - ENTER RESPECTIVE ITERATION (BASE & NODE) */
    if (rank == 0) {
        // mpi process of rank 0 in the MPI_COMM_WORLD will be the base station
        base_station(MPI_COMM_WORLD, new_comm, interval, nports);      
    } else {
        // other mpi process will be EV charging node
        charging_node(MPI_COMM_WORLD, new_comm, nrows, ncols, interval, nports, dims);     
    }

	MPI_Type_free( &AlertType );    /* Clean up the type */
    MPI_Finalize();
    return 0;
}

/**
 * @brief Base Station function handler
 * 
 * @param world_comm MPI_COMM_WORLD
 * @param comm new communicator for base station group
 * @param interval interval per iteration
 * @param nports number of ports
 * @return int 
 */
int base_station(MPI_Comm world_comm, MPI_Comm comm, int interval, int nports) {

    /* SECTION 0 - SETUP FOR LOGGING DATA */
    double simulation_start_time = MPI_Wtime();
    double simulation_end_time;
    double simulation_duration;

    char simulation_start_tstamp[20];
    char simulation_end_tstamp[20];

    char *temp = get_current_datetime();
    strcpy(simulation_start_tstamp, temp);
    
    int total_alert_detected = 0;
    int total_msg_exchanged = 0;
    double total_communication_time = 0;
    FILE *file = fopen("log_base.txt", "w");

    /* SECTION 1 - INITIAL SETUP FOR BASE STATION */
    int world_rank, world_size;
    int curr_iter = 0, pack_size, position;

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Pack_size(3, MPI_INT, MPI_COMM_WORLD, &pack_size);

    MPI_Request send_request[world_size-1];
    MPI_Status recv_status[world_size-1];
    MPI_Status send_status[world_size-1];

    int **charging_node_coords = NULL;
    int recv_cart_rank;
    char coord_recv_buffer[pack_size];
    charging_node_coords = (int**) malloc((world_size-1) * sizeof(int*));

    // Receive and store the coordinates of each charging node
    for (int i=1; i < world_size; i++) {
        position = 0;
        MPI_Recv(coord_recv_buffer, pack_size, MPI_PACKED, i, MSG_CART_DATA, world_comm, &recv_status[i-1]);
        MPI_Unpack(coord_recv_buffer, pack_size, &position, &recv_cart_rank, 1, MPI_INT, world_comm);
        
        charging_node_coords[recv_cart_rank] = (int*) malloc(2 * sizeof(int));
        MPI_Unpack(coord_recv_buffer, pack_size, &position, charging_node_coords[recv_cart_rank], 2, MPI_INT, world_comm);
    }

    pthread_t commThread;
    base_thread_data threadData;
    int *thread_shut_down;
    int node_termination_flag = 0;
    int num_msg_exchg_base_node[world_size-1];

    double *report_recv_time;
    int *report_node_arr_world;
    int *report_node_arr_cart;
    alert_data *report_detail_arr;
    int *num_report_recv;
    int curr_handling_idx;

    // Shared variables used between node and its ports
    report_recv_time = (double*) malloc((world_size-1) * sizeof(double));
    report_node_arr_world = (int*) malloc((world_size-1) * sizeof(int));
    report_node_arr_cart = (int*) malloc((world_size-1) * sizeof(int));;
    report_detail_arr = (alert_data*) malloc((world_size-1) * sizeof(*report_detail_arr));
    num_report_recv = (int*) malloc(sizeof(int));
    thread_shut_down = (int*) malloc(sizeof(int));
    *thread_shut_down = 0;

    /* SECTION 2 - WORKING ITERATION */
    while (curr_iter < MAX_ITERATION) {
        double start_time = MPI_Wtime();
        double end_time = 0;

        curr_handling_idx = 0;
        *num_report_recv = 0;

        // Setup the thread for communication with Charging Nodes
        if (curr_iter == 0) {
            threadData.report_recv_time = report_recv_time;
            threadData.report_node_arr_world = report_node_arr_world;
            threadData.report_node_arr_cart = report_node_arr_cart;
            threadData.report_detail_arr = report_detail_arr;
            threadData.num_report_recv = num_report_recv;
            threadData.shut_down = thread_shut_down;
            pthread_create(&commThread, NULL, ThreadFuncBaseRecv, &threadData);
        }

        // Reset the counting array
        for (int i=0; i < world_size-1; i++) {
            num_msg_exchg_base_node[i] = 0;
        }

        // INNER LOOP
        while (end_time - start_time < interval) {

            // If there exists alert/report that haven't been handled (handle it)
            if (curr_handling_idx != *num_report_recv) {
                alert_data curr_reporting_data = report_detail_arr[curr_handling_idx];
                int report_node_cart_rk = curr_reporting_data.my_node_rank;
                int report_node_world_rk = report_node_arr_world[curr_handling_idx];
                num_msg_exchg_base_node[report_node_world_rk-1] += 1;

                sleep(0.5);    // To provide time for report to be received (500 ms)

                int latest_report_idx = *num_report_recv;   // do this to avoid race condition with latest idx of report from base thread
                int free_nb_nodes[8];           // array of free nearby nodes (set to 8 as maximum of 8 possible free one) (CART RANK)
                int num_free_nb_nodes = 0;      // number of free nearby nodes
                int is_exists;
                int curr_nb_node;
                int direct_neigh_node;
                int neighbour_node_tracker_idx = 0;     // Keep track of the neighbour node of current report node that own the nearby node
                int num_valid_nb_nodes = 0;
                int valid_nb_nodes[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
                int is_new;

                /* Check if there is any free nearby nodes */
                for (int i=0; i < MAX_NEARBY_NEIGH; i++) {
                    is_exists = 0;
                    is_new = 1;
                    curr_nb_node = curr_reporting_data.nearby_nodes[i];
                    direct_neigh_node = curr_reporting_data.adjacent_nodes[neighbour_node_tracker_idx];

                    // No more valid nearby neighbours
                    if (curr_nb_node == -10) {
                        break;
                    } 

                    // If the nearby node valid (exist or not the reporting node itself)
                    else if (curr_nb_node != report_node_cart_rk && curr_nb_node != -2) {
                        for (int j=0; j < num_free_nb_nodes; j++) {
                            if (curr_nb_node == free_nb_nodes[j]) {
                                is_exists = 1;
                                break; 
                            }
                        }

                        if (!is_exists) {
                            // Check whether current nearby node is new
                            for (int i=0; i < num_valid_nb_nodes; i++) {
                                if (curr_nb_node == valid_nb_nodes[i]) {
                                    is_new = 0;
                                    break;
                                }
                            }

                            // It's new, record it
                            if (is_new) {
                                valid_nb_nodes[num_valid_nb_nodes] = curr_nb_node;
                                num_valid_nb_nodes += 1;
                            }

                            // Check whether the report received (both the nearby node itself and the neighbour node associated with it)
                            for (int j=0; j < latest_report_idx; j++) {
                                if (curr_nb_node == report_node_arr_cart[j] || direct_neigh_node == report_node_arr_cart[j]) {
                                    // if report received (move to next loop)
                                    break;
                                }  
                                else if (j == latest_report_idx-1) {  
                                    // if report not receive (up until then)     
                                    free_nb_nodes[num_free_nb_nodes] = curr_nb_node;
                                    num_free_nb_nodes++;
                                }
                            }
                        }
                    }
                    if ((i+1) % 4 == 0)
                        neighbour_node_tracker_idx += 1;
                }

                char resolution[256];
        
                if (num_free_nb_nodes > 0) {
                    // There is free nearby nodes
                    snprintf(resolution, sizeof(resolution), "Nearby Nodes Found! Please proceed to Node { ");
                    for (int i=0; i < num_free_nb_nodes-1; i++) {
                        int remaining_len = sizeof(resolution) - strlen(resolution);
                        snprintf(resolution + strlen(resolution), remaining_len, "%d,", free_nb_nodes[i]);
                    }
                    snprintf(resolution + strlen(resolution), sizeof(resolution) - strlen(resolution), "%d }\n", free_nb_nodes[num_free_nb_nodes-1]);
                    MPI_Isend(resolution, 256, MPI_CHAR, report_node_world_rk, MSG_RESOLUTION, world_comm, &send_request[0]);

                } else {
                    // No free nearby nodes
                    sprintf(resolution, "NO neighbour or nearby nodes available! Please wait for a while.\n");
                    MPI_Isend(resolution, 256, MPI_CHAR, report_node_world_rk, MSG_RESOLUTION, world_comm, &send_request[0]);
                }

                /* LOGGING ALERT DATA */
                char *alert_reported_time = curr_reporting_data.detected_timestamp;    
                char logged_time[20];
                temp = get_current_datetime();
                strcpy(logged_time, temp);
                free(temp);

                // number of adjacent node
                int num_adjacent_node = 0;
                for (int i=0; i < 4; i++) {
                    if (curr_reporting_data.adjacent_nodes[i] >= 0)
                        num_adjacent_node++;
                }

                int coord0 = charging_node_coords[report_node_cart_rk][0];  // row
                int coord1 = charging_node_coords[report_node_cart_rk][1];  // column
                int free_port = curr_reporting_data.num_free_port;          // number of free ports
                char *ip_add = curr_reporting_data.IPv4_add;                // ipv4 address

                // Header Row
                fprintf(file, "-------------------------------------------------------------------------------------------------------------------------------\n");
                fprintf(file, "Iteration: %d\n", curr_iter);
                fprintf(file, "%-35s%-25s\n", "Logged time :", alert_reported_time);
                fprintf(file, "%-35s%-25s\n", "Alert Reported time :", logged_time);
                fprintf(file, "Number of adjacent node: %d\n", num_adjacent_node);
                fprintf(file, "Number of nearby node: %d\n", num_valid_nb_nodes);
                fprintf(file, "Availability to be considered full: %d\n", MIN_PORT_AVAILABLE);

                fprintf(file, "\n");

                // Reporting Node 
                fprintf(file, "%-20s%-13s%-15s%-20s%-17s\n", "Reporting Node", "Coord", "Total Port", "Available Port", "IPv4");
                fprintf(file, "%-20d(%d, %d)%-7s%-15d%-20d%-17s\n", report_node_cart_rk, coord0, coord1, "", nports, free_port, ip_add);

                fprintf(file, "\n");

                // Adjacent Node 
                fprintf(file, "%-20s%-13s%-15s%-20s%-17s\n", "Adjacent Nodes", "Coord", "Total Port", "Available Port", "IPv4");
                for (int i=0; i < num_adjacent_node; i++) { 
                    int aj_node_rank = curr_reporting_data.adjacent_nodes[i];
                    coord0 = charging_node_coords[aj_node_rank][0];
                    coord1 = charging_node_coords[aj_node_rank][1];
                    free_port = curr_reporting_data.adjacent_node_port_usage[i];
                    ip_add = curr_reporting_data.neighbour_IPv4[i];
                    fprintf(file, "%-20d(%d, %d)%-7s%-15d%-20d%-17s\n", aj_node_rank, coord0, coord1, "", nports, free_port, ip_add);
                }

                fprintf(file, "\n");

                // Nearby Node
                fprintf(file, "%-20s%-13s\n", "Nearby Nodes", "Coord");
                for (int i=0; i < num_valid_nb_nodes; i++) { 
                    int nb_node_rank = valid_nb_nodes[i];
                    coord0 = charging_node_coords[nb_node_rank][0];
                    coord1 = charging_node_coords[nb_node_rank][1];
                    fprintf(file, "%-20d(%d, %d)%-7s\n", nb_node_rank, coord0, coord1, "");
                }

                fprintf(file, "\n");

                // Resolution
                fprintf(file, "Available station nearby (no report received in last 0.5 seconds): ");
                if (num_free_nb_nodes > 0) {
                    for (int i=0; i < num_free_nb_nodes-1; i++) {
                        fprintf(file, "%d, ", free_nb_nodes[i]);
                    }
                    fprintf(file, "%d\n", free_nb_nodes[num_free_nb_nodes-1]);
                } else {
                    fprintf(file, "None\n");
                }

                // Communication Time
                double time_offset = start_time - curr_reporting_data.node_start_time;
                double communication_time = report_recv_time[curr_handling_idx] - curr_reporting_data.sending_time - time_offset;
                total_communication_time += fabs(communication_time);
                fprintf(file, "Communication Time (seconds) : %lf\n", fabs(communication_time));

                // Total messages between base and reporting node (this iteration)
                int total_msg = num_msg_exchg_base_node[report_node_world_rk-1] + 1;
                fprintf(file, "Total Messages send between reporting node and base station: %d\n", total_msg);

                fprintf(file, "\n\n");

                curr_handling_idx++;
                total_msg_exchanged += (curr_reporting_data.num_msg_exchanged + total_msg);
                total_alert_detected++;
            }

            // Termination message for charging nodes
            if (curr_iter == MAX_ITERATION-2 && !node_termination_flag) {
                for (int i=1; i < world_size; i++) {
                    MPI_Isend(&curr_iter, 0, MPI_INT, i, MSG_EXIT, world_comm, &send_request[i-1]);
                }
                MPI_Waitall(world_size-1, send_request, send_status);
                node_termination_flag = 1;
            }

            end_time = MPI_Wtime();
            
        }

        // Termination message for the base SendRecv thread
        if (curr_iter == MAX_ITERATION-1) {
            *thread_shut_down = 1;
        }

        curr_iter += 1;
        MPI_Barrier(world_comm);
    }

    /* SECTION 3 - LOGGING SUMMARY DATA */
    temp = get_current_datetime();
    strcpy(simulation_end_tstamp, temp);

    simulation_end_time = MPI_Wtime();
    simulation_duration = simulation_end_time - simulation_start_time;

    fprintf(file, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    fprintf(file, "===============================================================================================================================\n");
    fprintf(file, "\n");

    fprintf(file, "SIMULATION SUMMARY\n\n");

    fprintf(file, "%-33s%-20s\n", "Simulation Start Timestamp:", simulation_start_tstamp);
    fprintf(file, "%-33s%-20s\n", "Simulation End Timestamp:", simulation_end_tstamp);
    fprintf(file, "%-33s%-20lf\n", "Simulation Duration (seconds):", simulation_duration);
    fprintf(file, "%-33s%-20d\n", "Total Iteration:", MAX_ITERATION);
    fprintf(file, "%-33s%-20d\n", "Number of Alert Detected:", total_alert_detected);
    fprintf(file, "%-33s%-20d\n", "Total Number of Messages:", total_msg_exchanged);
    fprintf(file, "%-33s%-20lf\n", "Total communication time:", total_communication_time);


    fprintf(file, "\n");
    fprintf(file, "===============================================================================================================================\n");
    fprintf(file, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");



    /* SECTION 4 - MEMORY CLEANING */
    for (int i=0; i < world_size-1; i++) {
        free(charging_node_coords[i]);
    }

    free(charging_node_coords);
    free(thread_shut_down);
    free(report_recv_time);
    free(report_node_arr_world);
    free(report_node_arr_cart);
    free(report_detail_arr);
    free(num_report_recv);

    return 0;
}


/**
 * @brief Charging Node Function Handler
 * 
 * @param world_comm MPI_COMM_WORLD
 * @param comm new communicator for the group of charging nodes
 * @param nrows number of rows
 * @param ncols number of cols
 * @param interval interval per iteration
 * @param nports number of ports
 * @param dims dimension
 * @return int 
 */
int charging_node(MPI_Comm world_comm, MPI_Comm comm, int nrows, int ncols, int interval, int nports, int *dims) {

    /* SECTION 0 - SETTING UP LOGGING DATA */
    char hostbuffer[256];
    char *IPbuffer;
    struct hostent *host_entry;

    gethostname(hostbuffer, sizeof(hostbuffer));     // To retrieve hostname
    host_entry = gethostbyname(hostbuffer);          // To retrieve host information
    IPbuffer = inet_ntoa(*((struct in_addr*)         // To convert an Internet network address into ASCII string
        host_entry->h_addr_list[0]));

    double total_comm_time_nodes = 0;       // total communication time between nodes 
    double total_comm_time_alert = 0;       // total communication time between node and station
    int total_not_full = 0;                 // number of times the node has available ports
    int total_full_with_neigh = 0;          // number of times the node has full ports but with available neighbour nodes
    int total_full_alert = 0;               // number of times the node has to alert the base station

    /* SECTION 1 - SETUP VIRTUAL CARTERSIAN TOPOLOGY */
    int ndims = 2, size, my_rank, world_size, my_cart_rank, reorder, ierr;
    int nbr_i_lo, nbr_i_hi;
    int nbr_j_lo, nbr_j_hi;
    int coord[ndims], wrap_around[ndims];
    MPI_Comm comm2D;
    int pack_size;
    int position = 0;
    
    MPI_Comm_size(world_comm, &world_size);    // size of the world communicator
    MPI_Comm_size(comm, &size);                // size of the new communicator (charging node)
    MPI_Comm_rank(comm, &my_rank);             // rank within new communicator (charging node)
    MPI_Dims_create(size, ndims, dims);        // create divison of processors in cartesian grid

    if (my_rank == 0) {
        printf("Root Rank: %d. Comm Size: %d. Grid Dimension = [%d x %d]\n", my_rank, size, dims[0], dims[1]);
    }

    /* create cartesian mapping */
    wrap_around[0]=wrap_around[1] = 0;     // periodic shift is false (both directions)
    reorder = 0;
    ierr = 0;

    MPI_Cart_create(comm, ndims, dims, wrap_around, reorder, &comm2D);  // make a new communicator on the topology info
    if (ierr != 0)
        printf("ERROR[%d] creating CART\n", ierr);

    /* find my coordinates in the cartesian communicator group */
    MPI_Cart_coords(comm2D, my_rank, ndims, coord);

    /* use my cartesian coordinates to find my rank in cartesian group*/
    MPI_Cart_rank(comm2D, coord, &my_cart_rank);

    /* get my neighbors; axis is coordinate dimension of shift */
    MPI_Cart_shift(comm2D, SHIFT_ROW, DISP, &nbr_i_lo, &nbr_i_hi);
    MPI_Cart_shift(comm2D, SHIFT_COL, DISP, &nbr_j_lo, &nbr_j_hi);

    MPI_Pack_size(3, MPI_INT, MPI_COMM_WORLD, &pack_size);
    char node_data_buf[pack_size];

    MPI_Pack(&my_cart_rank, 1, MPI_INT, node_data_buf, pack_size, &position, world_comm);
    MPI_Pack(coord, 2, MPI_INT, node_data_buf, pack_size, &position, world_comm);
    MPI_Send(node_data_buf, pack_size, MPI_PACKED, 0, MSG_CART_DATA, world_comm);


    /* SECTION 2 - INITIAL SETUP FOR THE CHARGING NODES */
    int neighbours[4] = {nbr_i_lo, nbr_i_hi, nbr_j_lo, nbr_j_hi};   // the neighbours of the current node
    MPI_Request send_request[4];
    MPI_Status send_status[4];
    MPI_Status status;
    MPI_Request request;
    MPI_Status status_world;

    // shared variables used between node and its ports
    int *curr_record_idx = (int*) malloc(sizeof(int));      // index for current position of shared array
    int *curr_iter = (int*) malloc(sizeof(int));            // current iteration number
    int *port_shut_down = (int*) malloc(sizeof(int));       // shut down indication for ports
    *curr_record_idx = 0;
    *curr_iter = 0;
    *port_shut_down = 0;

    pthread_t hThread[nports];
    node_thread_data threadData[nports];
    pthread_mutex_init(&shared_arr_Mutex, NULL);
    pthread_cond_init(&shared_arr_Cond, NULL);

    int is_full = 0;
    int neighbour_port_usage[4] = {-1, -1, -1, -1};     // neighbour node's port use condition (0 - available, 1 - full)
    int num_valid_neighbour = 0;
    int shut_down = 0;

    char fn_buf[256];
    FILE *pFile;
    sprintf(fn_buf, "log_%d.txt", my_cart_rank);
    pFile = fopen(fn_buf, "w");

    // Allocate memory for ports_record (NEED TO BE FREE)
    for (int i=0; i < MAX_SARRAY_SIZE; i++) {
        ports_record_arr[i] = (ports_record*) malloc(sizeof(ports_record));
        ports_record_arr[i]->total_free_port = 0;
    }

    /* SECTION 3 - WORKING ITERATION */
    while (true) {    

        if (*port_shut_down && shut_down) {
            // port already shut down and node is signalled to shutdown as well
            break;
        } else if (shut_down) {
            // port have not shut down (signal it)
            *port_shut_down = 1;
        }
        

        double start_time = MPI_Wtime();
        double end_time = 0;
        int num_msg_exchg_neigh = 0;
        is_full = 0;


        // Initialise the threads (ports)
        if (*curr_iter == 0) {
            for (int i=0; i < nports; i++) {
                threadData[i].node_rank = my_cart_rank;
                threadData[i].thread_rank = i;
                threadData[i].curr_record_idx = curr_record_idx;
                threadData[i].port_shut_down = port_shut_down;
                pthread_create(&hThread[i], NULL, ThreadFuncUpdatePort, &threadData[i]);
                sleep(0.8);
            }
        }

        pthread_mutex_lock(&shared_arr_Mutex);  

        // reset the occupied data (to accept new data)
        if (*curr_iter >= MAX_SARRAY_SIZE) {
            int reset_idx = *curr_iter % MAX_SARRAY_SIZE;
            free(ports_record_arr[reset_idx]->rec_timestamp);
            ports_record_arr[reset_idx]->total_free_port = 0;   
        }

        pthread_mutex_unlock(&shared_arr_Mutex);
        pthread_cond_broadcast(&shared_arr_Cond);

        char* curr_timestamp = get_current_datetime();
        ports_record_arr[*curr_record_idx]->rec_timestamp = strdup(curr_timestamp);
        free(curr_timestamp);

        sleep(0.8);
        MPI_Barrier(comm);

        // Calculate the size of packed data
        MPI_Pack_size(5, MPI_INT, MPI_COMM_WORLD, &pack_size);
        pack_size += 16 * sizeof(char);

        // Pack the data (the node's port usage & neighbour info)
        char send_buffer[pack_size];
        char recv_buffer[pack_size];
        position = 0;
        int num_ports_available = ports_record_arr[*curr_record_idx]->total_free_port;
        MPI_Pack(&num_ports_available, 1, MPI_INT, send_buffer, pack_size, &position, comm2D);
        MPI_Pack(&neighbours, 4, MPI_INT, send_buffer, pack_size, &position, comm2D);
        MPI_Pack(IPbuffer, 16, MPI_CHAR, send_buffer, pack_size, &position, comm2D);

        fprintf(pFile, "[ITERATION %d]\n", *curr_iter);
        fprintf(pFile, "Timestamp: %s\n", ports_record_arr[*curr_record_idx]->rec_timestamp);
        fprintf(pFile, "Available Ports: %d\n", num_ports_available);

        double requests_send_time = 0;
        double requests_recv_time = 0;
        double alert_send_time = 0;
        double alert_recv_time = 0;

        // Check if the ports is full / heavy utilised
        if (num_ports_available <= MIN_PORT_AVAILABLE) 
            is_full = 1;

        // Request for neighbour data (if full/heavy use)
        if (is_full) {
            requests_send_time = MPI_Wtime();
            num_valid_neighbour = 0;
            for (int i=0; i < MAX_NEIGHBOURS; i++) {
                send_request[i] = MPI_REQUEST_NULL; // reset request state
                int curr_neigh = neighbours[i];     // get current neighbour

                if (curr_neigh != -2) {
                    // if the neighbour exist, send to it
                    MPI_Isend(&is_full, 0, MPI_INT, curr_neigh, MSG_REQUEST_USAGE, comm2D, &send_request[i]);
                    num_valid_neighbour++;
                    num_msg_exchg_neigh++;
                }
            }
            MPI_Waitall(num_valid_neighbour, send_request, send_status);
        } else {
            total_not_full += 1;
            fprintf(pFile, "Condition: NOT FULL\n\n");
        }

        int num_neigh_data_recv = 0;
        int nearby_nodes[16];
        int curr_recv_pos = 0;
        int neigh_data_recv_order[MAX_NEIGHBOURS] = {-2, -2, -2, -2};
        int free_neighbour[MAX_NEIGHBOURS] = {-2, -2, -2, -2};
        char neighbour_ipv4[MAX_NEIGHBOURS][16];

        // initialise 
        for (int i=0; i < MAX_NEARBY_NEIGH; i++) {
            nearby_nodes[i] = -10;
        }

        // inner loop (send, recv and do operation)
        while (end_time - start_time < interval) {
            int flag2D = 0;
            int flagWorld = 0;
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm2D, &flag2D, &status);

            MPI_Iprobe(0, MPI_ANY_TAG, world_comm, &flagWorld, &status_world);
            
            if (flag2D && status.MPI_SOURCE != my_cart_rank) {
                // If there is incoming message from neighbouring nodes (charging nodes) and it's not the node itself

                switch (status.MPI_TAG)
                {
                    case MSG_REQUEST_USAGE:  /* Send current node's data to requester */
                        MPI_Recv(&is_full, 0, MPI_INT, MPI_ANY_SOURCE, MSG_REQUEST_USAGE, comm2D, &status);
                        MPI_Isend(send_buffer, pack_size, MPI_PACKED, status.MPI_SOURCE, MSG_SENT_USAGE, comm2D, &request);
                        break;

                    case MSG_SENT_USAGE:    /* Handling received neighbour data */
                        position = 0;
                        MPI_Recv(recv_buffer, pack_size, MPI_PACKED, MPI_ANY_SOURCE, MSG_SENT_USAGE, comm2D, &status);
                        MPI_Unpack(recv_buffer, pack_size, &position, &neighbour_port_usage[num_neigh_data_recv], 1, MPI_INT, comm2D);
                        MPI_Unpack(recv_buffer, pack_size, &position, &nearby_nodes[curr_recv_pos], 4, MPI_INT, comm2D);
                        MPI_Unpack(recv_buffer, pack_size, &position, &neighbour_ipv4[num_neigh_data_recv], 16, MPI_CHAR, comm2D);

                        neigh_data_recv_order[num_neigh_data_recv] = status.MPI_SOURCE;
                        curr_recv_pos += 4;
                        num_neigh_data_recv += 1;
                        num_msg_exchg_neigh += 1;

                        if (num_neigh_data_recv == num_valid_neighbour) {
                            // When receive all the requested neighbour data
                            requests_recv_time = MPI_Wtime();
                            total_comm_time_nodes += fabs(requests_recv_time - requests_send_time);         // accumulate the communication time
                            int num_free_node = 0;
                            for (int i=0; i < MAX_NEIGHBOURS; i++) {
                                if (neigh_data_recv_order[i] == -2) 
                                    break;

                                if (neighbour_port_usage[i] > MIN_PORT_AVAILABLE) {
                                    free_neighbour[num_free_node] = neigh_data_recv_order[i];
                                    num_free_node += 1;
                                } 
                            }
                            // Has FREE neighbour
                            if (num_free_node > 0) {
                                total_full_with_neigh += 1;
                                fprintf(pFile, "Condition: FULL with Available Neighbour Nodes\n");
                                fprintf(pFile, "Solution: Proceed to Node { ");

                                for (int i=0; i < num_free_node-1; i++) 
                                    fprintf(pFile, "%d,", free_neighbour[i]);

                                fprintf(pFile, "%d }\n\n", free_neighbour[num_free_node-1]);


                            } else {
                                // No free neighbour (alert and seek for resolution from base station)
                                total_full_alert += 1;
                                alert_data alert;
                                char *temp_timestamp = get_current_datetime();

                                alert_send_time = MPI_Wtime();
                                
                                /* Initialise all details needed for the alert */
                                strcpy(alert.detected_timestamp, temp_timestamp);
                                alert.num_msg_exchanged = num_msg_exchg_neigh;
                                alert.sending_time = alert_send_time;
                                alert.node_start_time = start_time;
                                alert.num_free_port = num_ports_available;
                                alert.my_node_rank = my_cart_rank;
                                strcpy(alert.IPv4_add, IPbuffer);
                                
                                for (int i=0; i < MAX_NEIGHBOURS; i++) {
                                    alert.adjacent_nodes[i] = neigh_data_recv_order[i];
                                    alert.adjacent_node_port_usage[i] = neighbour_port_usage[i];
                                    strcpy(alert.neighbour_IPv4[i], neighbour_ipv4[i]);
                                }

                                for (int i=0; i < MAX_NEARBY_NEIGH; i++) {
                                    alert.nearby_nodes[i] = nearby_nodes[i];
                                }

                                // Send the alert to the base station
                                MPI_Isend(&alert, 1, AlertType, 0, MSG_ALERT, world_comm, &request);
                                free(temp_timestamp);

                                fprintf(pFile, "Condition: ALERT BASE STATION\n");
                                MPI_Wait(&request, &status);
                            }
                        }
                    break;
                }

            } else if (flagWorld) {
                // If there is incoming message from base station

                switch (status_world.MPI_TAG)
                {
                    case MSG_EXIT:      /* End the Charging Node */
                        MPI_Recv(&is_full, 0, MPI_INT, 0, MSG_EXIT, world_comm, &status_world);
                        shut_down = 1;      // ready to shut down in the 2 loop after 
                        break;

                    case MSG_RESOLUTION:        /* Solution for full ports */ 
                        char resolution[256];
                        MPI_Recv(resolution, 256, MPI_CHAR, 0, MSG_RESOLUTION, world_comm, &status_world);
                        alert_recv_time = MPI_Wtime();
                        total_comm_time_alert += fabs(alert_recv_time - alert_send_time);
                        fprintf(pFile, "Solution: %s\n", resolution);
                        break;
                }
            }
            end_time = MPI_Wtime();
        }

        *curr_record_idx = ((*curr_record_idx + 1) % MAX_SARRAY_SIZE);
        *curr_iter += 1;

        MPI_Barrier(world_comm);    // Wait all nodes and base station reach the common point of end of iteration 
    }

    // Join
    for (int i = 0; i < nports; i++) {
        pthread_join(hThread[i], NULL);
    }

    /* SECTION 4 - LOGGING */
    fprintf(pFile, "\n");
    fprintf(pFile, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
    fprintf(pFile, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
    fprintf(pFile, "\n");

    fprintf(pFile, "SIMULATION SUMMARY\n\n");

    fprintf(pFile, "%-65s%-20d\n", "Number of Iteration with available ports: ", total_not_full);
    fprintf(pFile, "%-65s%-20d\n", "Number of Iteration with full ports (available neighbours): ", total_full_with_neigh);
    fprintf(pFile, "%-65s%-20d\n", "Number of Iteration with full ports (alert): ", total_full_alert);
    fprintf(pFile, "%-65s%-20lf\n", "Average communication time with neighbour nodes (seconds): ", total_comm_time_nodes/total_full_with_neigh);
    fprintf(pFile, "%-65s%-20lf\n", "Average communication time with base staiton (seconds): ", total_comm_time_alert/total_full_alert);

    fprintf(pFile, "\n");
    fprintf(pFile, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
    fprintf(pFile, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");

    /* SECTION 5 - MEMORY CLEANING */
    for (int i=0; i < MAX_SARRAY_SIZE; i++) {
        free(ports_record_arr[i]->rec_timestamp);
        free(ports_record_arr[i]);
    }

    free(curr_record_idx);
    free(curr_iter);

    pthread_mutex_destroy(&shared_arr_Mutex);
    pthread_cond_destroy(&shared_arr_Cond);

    return 0;
}


/**
 * @brief Get the current datetime in string format
 * 
 * @return char* current datetime 
 */
char* get_current_datetime() {
    char datetime_format[] = "%d/%m/%Y %H:%M:%S";  // format for the timestamp
    time_t rawtime;
    time(&rawtime);     // get the current time in seconds
    rawtime += 28800;   // add 8 hours (UTC -> MYT)

    struct tm *curr_dt;
    char* curr_dt_str = (char*) malloc(20 * sizeof(char));  // string representation of current datetime
    curr_dt = localtime( &rawtime );    // convert the time to tm strcuture format 

    strftime(curr_dt_str, 20, datetime_format, curr_dt);   // Convert the time into the defined format (string)

    return curr_dt_str;
}

/**
 * @brief Randomly update availability
 * 
 * @param th_rank thread rank
 * @param nd_rank node rank
 * @return int 1 or 0
 */
int rand_update_availabilty(int th_rank, int nd_rank) {
    
    srand(time(NULL) + th_rank*nd_rank);      // Seed the random number generator with current time, thread rank and node rank
    int is_available = rand() % 2;            // Randomly update availability of port ( 0 - in use, 1 - free)

    return is_available;
}


/**
 * @brief Thread function for port to update the availability
 * 
 * @param pArg thread data needed 
 * @return void* 
 */
void *ThreadFuncUpdatePort(void *pArg) {
    node_thread_data *p_td = (node_thread_data*) pArg;
    int th_rank = p_td->thread_rank;
    int nd_rank = p_td->node_rank;
    int *cr_record_idx = p_td->curr_record_idx;   
    int *pt_shut_down = p_td->port_shut_down; 

    while (true) {
        // Randomly update availability
        int curr_availability = rand_update_availabilty(th_rank, nd_rank);
        

        pthread_mutex_lock(&shared_arr_Mutex);
        pthread_cond_wait(&shared_arr_Cond, &shared_arr_Mutex);

        ports_record_arr[*cr_record_idx]->total_free_port += curr_availability;
        pthread_mutex_unlock(&shared_arr_Mutex);

        // Shut down the port or thread
        if (*pt_shut_down) {
            return NULL;
        }
    }

}

/**
 * @brief Thread function used for receiving the alert from charging nodes
 * 
 * @param pArg thread data needed (shared variables purposes)
 * @return void* 
 */
void *ThreadFuncBaseRecv(void *pArg) {
    base_thread_data *b_td = (base_thread_data*) pArg;
    double *rp_recv_time = b_td->report_recv_time;              // time where the report/alert is received by base 
    int *rp_node_arr_world = b_td->report_node_arr_world;       // rank of the node that make report (WORLD RANK)
    int *rp_node_arr_cart = b_td->report_node_arr_cart;         // rank of the node that make report (CART RANK)
    alert_data *rp_detail_arr = b_td->report_detail_arr;        // report detail of the alert
    int *num_rp_recv = b_td->num_report_recv;                   // number of report received in this iteration (up until now)
    int *th_shut_down = b_td->shut_down;                        // shut down flag

    alert_data alert_recv;
    MPI_Status status;

    while (true) {
        int flag_alert = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, MSG_ALERT, MPI_COMM_WORLD, &flag_alert, &status);

        if (flag_alert) {
            // There is a alert

            MPI_Recv(&alert_recv, 1, AlertType, MPI_ANY_SOURCE, MSG_ALERT, MPI_COMM_WORLD, &status); // Receive it
            
            /* Update all the details related to the alert */
            double recv_time = MPI_Wtime();
            rp_recv_time[*num_rp_recv] = recv_time;
            rp_node_arr_world[*num_rp_recv] = status.MPI_SOURCE; 
            rp_node_arr_cart[*num_rp_recv] = alert_recv.my_node_rank;   
            rp_detail_arr[*num_rp_recv] = alert_recv;    
            *num_rp_recv += 1;      
        }

        // Handling the shutdown condition 
        if (*th_shut_down) {
            return NULL;
        }
    }
}
