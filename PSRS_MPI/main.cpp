# include<mpi.h>
# include<bits/stdc++.h>

# include"../Utils.cpp"
# include"../Log.cpp"
using namespace std;

// 进程内部排序的时候使用qsort，需提供cmp函数，从小到大
int cmp(const void *a, const void *b) {
    return (*(int *)a - *(int *)b);
}

void swap(int& m, int& n){
    int temp = m;
    m = n;
    n = temp;
}

// 串行算法：标准库
void Sort_seq(vector<int>& nums, int begin, int end){
    sort(nums.begin() + begin, nums.begin() + end);
}

void Logdebug(int nums[], int n){
    for(int i = 0; i < n; ++i){
        cout << nums[i] << " ";
    }
    cout << endl;
}

void merge_two(int *array, int begin1, int end1, int begin2, int end2){
    int iter1 = begin1;
    int iter2 = begin2;
    int merged_iter = 0;
    int merged[end1 - begin1 + end2 - begin2];

    while(iter1 != end1 && iter2 != end2){
        if(array[iter1] > array[iter2]){
            merged[merged_iter] = array[iter2];
            merged_iter += 1;
            iter2 += 1;
        }
        else{
            merged[merged_iter] = array[iter1];
            merged_iter += 1;
            iter1 += 1;            
        }
    }
    while(iter1 != end1){
        merged[merged_iter] = array[iter1];
        iter1 += 1;
        merged_iter += 1;
    }
    while(iter2 != end2){
        merged[merged_iter] = array[iter2];
        merged_iter += 1;
        iter2 += 1;   
    }

    // 归回原位
    for(int i = 0; i < merged_iter; ++i){
        array[begin1 + i] = merged[i];
    }
    return;
}


void merge_sort(int *new_distbuf, int *recv_disp, int *recv_count, int n){
    if(n == 1)
        return;

    for(int i = 1; i < n; ++i){
        merge_two(new_distbuf, recv_disp[0], recv_disp[i], recv_disp[i], recv_disp[i] + recv_count[i]);
    }
}

// 正则采样排序
void PSRS(vector<int>& nums, int argc,char* argv[]){
    int proc_num; 
    int proc_id;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    // cout << "This is process " << proc_id << " of " << proc_num << endl;

    // 1. 均匀划分

    int block_size = nums.size() / proc_num;
    int block_size_pro = block_size + nums.size() % proc_num;

    // 计算各块的数据偏移displs与数据量sendcounts
    // 将无法整除的余数堆到最后一个进程，该块数据量为block_size_pro
    int sendcounts[proc_num];
    int displs[proc_num];
    for(int i = 0; i < proc_num; ++i){
        displs[i] = i * block_size;
        sendcounts[i] = block_size;
    }
    sendcounts[proc_num - 1] = block_size_pro;

    // 数据接收缓冲区
    int distbuf[block_size_pro];
    int* rootbuf = &nums[0];
    MPI_Scatterv(rootbuf, sendcounts, displs, MPI_INT, distbuf, block_size_pro, MPI_INT, 0, MPI_COMM_WORLD);
    // cout << "process "<< proc_id << " uniform divided" << endl;
    
    // 2. 局部排序
    if(proc_id == proc_num - 1)
        qsort(distbuf, block_size_pro, sizeof(MPI_INT), cmp);
    else
        qsort(distbuf, block_size, sizeof(MPI_INT), cmp);
    // cout << "process "<< proc_id << " locally sorted: ";
    // Logdebug(distbuf, block_size_pro);

    // 3. 选取样本
    // p个进程每个进程选取p个元素，采样数目等于进程数目
    // 即使是block_size_pro也不所谓，反正尾巴mod部分就当看不见
    int samples[proc_num * proc_num];
    int sample[proc_num];
    for(int i = 0; i < proc_num; ++i){
        sample[i] = distbuf[i * block_size / proc_num];
    }
    MPI_Gather(sample, proc_num, MPI_INT, samples, proc_num, MPI_INT, 0, MPI_COMM_WORLD);
    // cout << "process "<< proc_id << " sampling" << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    // 4. 样本排序(p^2个)
    if(proc_id == 0){
        qsort(samples, proc_num * proc_num, sizeof(MPI_INT), cmp);
        // cout << "main process" << " sort sample: ";
        //Logdebug(samples, 4);
    }

    // 5. 选择主元（p-1个）
    int primes[proc_num - 1];
    if(proc_id == 0){
        for(int i = 0; i < proc_num - 1; ++i){
            primes[i] = samples[(i + 1) * proc_num];
        }
    }
    // cout << "process "<< proc_id << " prime chosen" << endl;

    // 6. 主元划分
    MPI_Bcast(primes, proc_num - 1, MPI_INT, 0, MPI_COMM_WORLD);
    // cout << "prime chosen: " << primes[0] << endl;
    // 按照主元重塑各分段
    vector<vector<int>> reshapedNums(proc_num, vector<int>());

    // 需要区分pro与否
    if(proc_id == proc_num - 1){
        // 外层循环遍历各主元，每次构建一行
        int j = 0;
        for(int i = 0; i < proc_num - 1; ++i){
            int partition_prime = primes[i];
            for(;j < block_size_pro && distbuf[j] < partition_prime; ++j){
                reshapedNums[i].push_back(distbuf[j]);
            }
        }
        // 最后一行
        reshapedNums.back() = vector<int>(block_size_pro - j);
    }
    else{
        // 外层循环遍历各主元，每次构建一行
        int j = 0;
        for(int i = 0; i < proc_num - 1; ++i){
            int partition_prime = primes[i];
            // 万一一开始就除不尽那不是完蛋？真实数目没有block_size那么多
            for(;j < block_size && distbuf[j] < partition_prime; ++j){
                reshapedNums[i].push_back(distbuf[j]);
            }
        }        
        // 最后一行
        reshapedNums.back() = vector<int>(block_size - j);        
    }

    // 7. 全局交换
    // send_count[]存储各分段的大小, send_disp[]存储偏移量
    int send_count[proc_num];
    int send_disp[proc_num];
    send_disp[0] = 0;
    for(int i = 0; i < reshapedNums.size(); ++i){
        send_count[i] = reshapedNums[i].size();
        if(i < reshapedNums.size() - 1)
            send_disp[i + 1] = send_disp[i] + send_count[i];
    }
    // cout << "process "<< proc_id << " prime redivided, length is: ";
    //Logdebug(send_count, 2);

    // 在数据交换之前，还需进行数据长度的全局交换，计算recv_size与recv_count
    int recv_count[proc_num];
    int recv_disp[proc_num]; recv_disp[0] = 0;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MPI_COMM_WORLD);
    // recv_size累加收到数据的总数
    int recv_size = recv_count[0];

    for(int i = 1; i < proc_num; ++i){
        recv_disp[i] = recv_count[i - 1] + recv_disp[i - 1];
        recv_size += recv_count[i];
    }

    // 正式交换数据，使用新的buffer来接收
    int new_distbuf[recv_size];
    MPI_Alltoallv(distbuf, send_count, send_disp, MPI_INT, 
        new_distbuf, recv_count, recv_disp, MPI_INT, 
        MPI_COMM_WORLD
    );
    // cout << "process "<< proc_id << " globally exchanged" << endl;
    
    
    // 8. 归并排序
    merge_sort(new_distbuf, recv_disp, recv_count, proc_num);
    // cout << "process "<< proc_id << " resort: ";
    //Logdebug(new_distbuf, recv_size);

    MPI_Barrier(MPI_COMM_WORLD);

    // 9. 收集结果
    // 各处理器中的元素总数汇报一下，复用recv_count数组
    MPI_Gather(&recv_size, 1, MPI_INT, recv_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // cout << "main process "<< proc_id << " has nums size: ";
    //Logdebug(recv_count, 2);

    recv_disp[0] = 0;
    for(int i = 1; i < proc_num; ++i){
        recv_disp[i] = recv_disp[i - 1] + recv_count[i - 1];
    }

    // cout << "Main process receive nums at: ";
    // Logdebug(recv_disp, 2);
    // 将各处理器的结果汇总到rootbuf
    MPI_Gatherv(new_distbuf, recv_size, MPI_INT, rootbuf, recv_count, recv_disp, MPI_INT, 0, MPI_COMM_WORLD);
    // cout << "process "<< proc_id << " result merge" << endl;

    // 10. 释放资源
    MPI_Finalize();
}


void Main(int argc, char* argv[], string batch_size, Log* log, int proc_num){
    // 读数据
    vector<int> nums1;
    vector<int> nums2;
    string filename = "dataset/Data_" + batch_size + ".txt";
    readData(filename, nums1);
    readData(filename, nums2);
    int n = nums1.size();

    // 串行排序计时
    auto start = system_clock::now();
    Sort_seq(nums1, 0, n);
    auto end = system_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    double seq_time = double(duration.count() * microseconds::period::num) / milliseconds::period::den;

    // PSRS计时
    start = system_clock::now();
    PSRS(nums2, argc, argv);
    end = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    double PSRS_time = double(duration.count() * microseconds::period::num) / milliseconds::period::den;

    // 正确性检验，差分测试
    assert(CorrectnessTest(nums1, nums2));

    // 写入日志
    log->write(batch_size, seq_time, PSRS_time);
}


int main(int argc,char* argv[]){

    int proc_num = atoi(argv[2]);
    // 变量一：线程数，由程序执行时指定
    string filename = "Results/Result_T" + to_string(proc_num) + ".log";
    Log* log = new Log(filename);

    // 变量二：数据集大小
    // Main(argc, argv, "1K", log, proc_num);
    // Main(argc, argv, "5K", log, proc_num);
    // Main(argc, argv, "10K", log, proc_num);
    Main(argc, argv, "100K", log, proc_num);

    // 关闭文件以使得修改可见
    log->~Log();        

    return 0;
}
