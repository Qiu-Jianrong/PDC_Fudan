# include "../Log.cpp"
# include "../Utils.cpp"
# include<mpi.h>
# include<bits/stdc++.h>

using namespace std;

int dims[2], periods[2];
MPI_Comm comm_2d;// 二维环绕结构

// 普通矩阵相乘，可以考虑ssthresh
void Mul_seq(const vector<vector<float>>& A, const vector<vector<float>>& B, vector<vector<float>>& C, int& n){

    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            C[i][j] = 0;
            for(int k = 0; k < n; ++k){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void Mul_cannon(float* A, float* B, float* C, int n){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            C[i * n + j] = 0;
            for(int k = 0; k < n; ++k){
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
}


// 简单分块乘法
void Simple_block_para(int argc, char* argv[], vector<vector<float>>& A, vector<vector<float>>& B, vector<vector<float>>& C, int& n){
    MPI_Init(&argc, &argv);
    int proc_num;
    int proc_id;
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    int proc_num_inline = sqrt(proc_num);
    int block_size = n / proc_num_inline;

    
    float local_A[block_size * block_size];
    float local_B[block_size * block_size];
    int x_coor = proc_id / proc_num_inline;
    int y_coor = proc_id % proc_num_inline;
    // (i, j)进程得到[i * block_size ~ (i + 1) * block_size]*[j * block_size ~ (j + 1)* block_size]
    for(int i = 0; i < block_size; ++i){
        for(int j = 0; j < block_size; ++j){
            local_A[i * block_size + j] = A[x_coor * block_size + i][y_coor * block_size + j];
            local_B[i * block_size + j] = B[x_coor * block_size + i][y_coor * block_size + j];
        }
    }

    // 块(i, j)需要接收所有的A[i][k]和B[k][j]
    // 以正数表示A矩阵，负数表示B矩阵
    // 不必！每个块收到的源头都不一样
    for(int k = 0; k < proc_num_inline; ++k){
        if(proc_id != x_coor * proc_num_inline + k)
            MPI_Send(local_A, block_size * block_size, MPI_FLOAT, 
                x_coor * proc_num_inline + k, proc_id, MPI_COMM_WORLD);
        if(proc_id != k * proc_num_inline + y_coor)
            MPI_Send(local_B, block_size * block_size, MPI_FLOAT, 
                k * proc_num_inline + y_coor, proc_id, MPI_COMM_WORLD);
    }
    float recv_A[block_size * block_size * proc_num_inline];
    float recv_B[block_size * block_size * proc_num_inline];  
    float local_C[block_size * block_size];  
    MPI_Status status;
    for(int k = 0; k < proc_num_inline; ++k){
        if(proc_id != x_coor * proc_num_inline + k)
            MPI_Recv(recv_A + (k * block_size * block_size), block_size * block_size, MPI_FLOAT, 
                x_coor * proc_num_inline + k, x_coor * proc_num_inline + k, MPI_COMM_WORLD, &status);
        else
            memcpy(recv_A + (k * block_size * block_size), local_A, block_size * block_size * sizeof(float));
        if(proc_id != k * proc_num_inline + y_coor)
            MPI_Recv(recv_B + (k * block_size * block_size), block_size * block_size, MPI_FLOAT, 
                k * proc_num_inline + y_coor, k * proc_num_inline + y_coor, MPI_COMM_WORLD, &status);
        else
            memcpy(recv_B + (k * block_size * block_size), local_B, block_size * block_size * sizeof(float));
        Mul_cannon(recv_A + (k * block_size * block_size), 
            recv_B + (k * block_size * block_size), local_C, block_size);
    }
    
    // 收集结果
    for(int i = 0; i < block_size; ++i){
        for(int j = 0; j < block_size; ++j){
            C[x_coor * block_size + i][y_coor * block_size + j] = local_C[i * block_size + j];
        }
    }

    MPI_Finalize();
}

// Cannon矩阵相乘
void Cannon(int argc, char* argv[], vector<vector<float>>& A, vector<vector<float>>& B, vector<vector<float>>& C, int& n){
    int proc_num; 
    int proc_id;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    int proc_num_inline = sqrt(proc_num);
    // 每个矩阵块的行数大小，假设均可整除
    int block_size = n / proc_num_inline;


    // 使用2D环绕连接，适合循环移位！
    dims[0] = dims[1] = proc_num_inline;
    periods[0] = periods[1] = 0;// 0表示二维网孔结构
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_2d);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

    // 获取当前处理器在二维循环结构中的坐标，准备好发送与接收缓冲区并初始化
    int my_coor[2];
    MPI_Cart_coords(comm_2d, proc_id, 2, my_coor);
    float local_A[block_size * block_size];
    float local_B[block_size * block_size];
    // (i, j)进程得到[i * block_size ~ (i + 1) * block_size]*[j * block_size ~ (j + 1)* block_size]
    for(int i = 0; i < block_size; ++i){
        for(int j = 0; j < block_size; ++j){
            local_A[i * block_size + j] = A[my_coor[0] * block_size + i][my_coor[1] * block_size + j];
            local_B[i * block_size + j] = B[my_coor[0] * block_size + i][my_coor[1] * block_size + j];
        }
    }

    // 1. 第一次重排列，进程(i, j)发送A到进程
    // (i, (j + proc_num_inline - i) % proc_num_inline)
    // 进程(i, j)发送B到((i + proc_num_inline - j) % proc_num_inline, j)
    // tag就采用发送方的rank
    int left_pos = my_coor[0] * proc_num_inline + 
        ((my_coor[1] + proc_num_inline - my_coor[0]) % proc_num_inline);
    int up_pos = ((my_coor[0] + proc_num_inline - my_coor[1]) % proc_num_inline) 
        * proc_num_inline + my_coor[1];
    // int source_pos = my_coor[0] * proc_num_inline + my_coor[1];
    MPI_Status status;
    MPI_Sendrecv_replace(local_A, block_size * block_size, MPI_FLOAT, 
        left_pos, proc_id, proc_id, proc_id, comm_2d, &status);
    MPI_Sendrecv_replace(local_B, block_size * block_size, MPI_FLOAT, 
        up_pos, proc_id, proc_id, proc_id, comm_2d, &status);

        
    float local_C[block_size * block_size];
    // 2. 所有处理器并行乘加
    Mul_cannon(local_A, local_B, local_C, block_size);

    left_pos = my_coor[0] * proc_num_inline + 
        ((my_coor[1] + proc_num_inline - 1) % proc_num_inline);
    up_pos = ((my_coor[0] + proc_num_inline - 1) % proc_num_inline) 
        * proc_num_inline + my_coor[1];
    // 循环proc_num_inline - 1次
    for(int i = 0; i < proc_num_inline - 1; ++i){
        // 3. 向上/向左一格，所有处理器并行乘加
        MPI_Sendrecv_replace(local_A, block_size * block_size, MPI_FLOAT, 
            left_pos, proc_id, proc_id, proc_id, comm_2d, &status);
        MPI_Sendrecv_replace(local_B, block_size * block_size, MPI_FLOAT, 
            up_pos, proc_id, proc_id, proc_id, comm_2d, &status);
        Mul_cannon(local_A, local_B, local_C, block_size);
    }
    
    // 4. 收集结果并变换形式
    for(int i = 0; i < block_size; ++i){
        for(int j = 0; j < block_size; ++j){
            C[my_coor[0] * block_size + i][my_coor[1] * block_size + j] = local_C[i * block_size + j];
        }
    }

    MPI_Finalize();
}

void Main(int argc, char* argv[], int n, Log* log){
    // 读数据
    vector<vector<float>> A(n, vector<float>(n, 0));
    vector<vector<float>> B(n, vector<float>(n, 0));
    string filename = "dataset/Data_" + to_string(n) + ".txt";
    readMatrix(filename, A, n);
    readMatrix(filename, B, n);
    vector<vector<float>> C1(n, vector<float>(n, 0));
    vector<vector<float>> C2(n, vector<float>(n, 0));

    // 串行计时
    auto start = system_clock::now();
    Mul_seq(A, B, C1, n);
    auto end = system_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    double seq_time = double(duration.count() * microseconds::period::num) / milliseconds::period::den;

    // 简单分块并行计时
    // start = system_clock::now();
    // Simple_block_para(argc, argv, A, B, C2, n);
    // end = system_clock::now();
    // duration = duration_cast<microseconds>(end - start);
    // double simple_block_time = double(duration.count() * microseconds::period::num) / milliseconds::period::den;

    // cannon计时
    start = system_clock::now();
    Cannon(argc, argv, A, B, C2, n);
    end = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    double cannon_time = double(duration.count() * microseconds::period::num) / milliseconds::period::den;


    // 正确性检验，差分测试
    assert(MatrixCorrectness(C1, C2, n));

    // 写入日志
    log->writeMatrix(n, seq_time, cannon_time);
    // log->writeMatrix(n, seq_time, simple_block_time);
}


int main(int argc, char* argv[]){
    // int batch_size = 32;
    // int batch_size = 64;
    // int batch_size = 128;
    int batch_size = 256;
    string filename = "Results/Result_T" + to_string(stoi(argv[2])) + ".log";
    Log* log = new Log(filename);
    Main(argc, argv, batch_size, log);
    log->~Log();
}