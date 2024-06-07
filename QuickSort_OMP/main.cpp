# include<bits/stdc++.h>
# include<omp.h>

# include"../Utils.cpp"
#include"../Log.cpp"
using namespace std;

void swap(int& m, int& n){
    int temp = m;
    m = n;
    n = temp;
}

void QS_seq(vector<int>& nums, int begin, int end){
    if(begin == end)
        return;

    // 随机化
    swap(nums[begin], nums[(end + begin) / 2]);

    int pivot = nums[begin];
    int pivot_pos = begin;
    for(int i = begin + 1; i < end; ++i){
        if(nums[i] < pivot){
            swap(nums[pivot_pos + 1], nums[i]);
            swap(nums[pivot_pos], nums[pivot_pos + 1]);
            pivot_pos += 1;
        }
    }

    QS_seq(nums, begin, pivot_pos);
    QS_seq(nums, pivot_pos + 1, end);
}

// 并行快排
void QS_para(vector<int>& nums, int begin, int end, int cutoff){
    if(begin == end)
        return;
    

    // 随机化
    swap(nums[begin], nums[(end + begin) / 2]);

    // 分割数组
    int pivot = nums[begin];
    int pivot_pos = begin;
    for(int i = begin + 1; i < end; ++i){
        if(nums[i] < pivot){
            swap(nums[pivot_pos + 1], nums[i]);
            swap(nums[pivot_pos], nums[pivot_pos + 1]);
            pivot_pos += 1;
        }
    }

    // 递归排序子数组，如果数组长度足够小，改用串行
    if(end - begin <= cutoff){
        QS_seq(nums, begin, pivot_pos);
        QS_seq(nums, pivot_pos + 1, end);
        return;
    }

    #pragma omp task shared(nums)
    QS_para(nums, begin, pivot_pos, cutoff);

    #pragma omp task shared(nums)
    QS_para(nums, pivot_pos + 1, end, cutoff);        
    
    #pragma omp taskwait
}

// 测试数据并写入日志
void Main(const string& filename, const string& batch, Log* log, int threads_num, void (*Sort_seq)(vector<int>&, int, int), void (*Sort_para)(vector<int>&, int, int, int)){
    vector<int> nums1, nums2;
    // string filename = "dataset/Data_" + batch + ".txt";
    readData(filename, nums1);
    readData(filename, nums2);

    int n = nums1.size();
    int cutoff = 100;

    // 串行快排计时
    auto start = system_clock::now();
    // _sleep(1000);
    Sort_seq(nums1, 0, n);
    auto end = system_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    double seq_time = double(duration.count() * microseconds::period::num) / milliseconds::period::den;

    // 并行快排计时
    start = system_clock::now();

    #pragma omp parallel num_threads(threads_num)
    {
        #pragma omp single
        Sort_para(nums2, 0, n, cutoff);
    }
    
    end = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    double para_time = double(duration.count() * microseconds::period::num) / milliseconds::period::den;

    // 正确性检验
    assert(CorrectnessTest(nums1, nums2));

    // 写入日志
    log->write(batch, seq_time, para_time);
}

int main(){
    srand(time(NULL));

    // 变量一：线程数
    for(int i = 2; i < 17; i += 2){
        string filename = "Results/Result_T" + to_string(i) + ".log";
        Log* log = new Log(filename);
        // 变量二：数据集大小
        // 因为是工作目录所以不用加前缀
        Main("dataset/Data_1K.txt", "1K", log, i, QS_seq, QS_para);
        Main("dataset/Data_5K.txt", "5K", log, i, QS_seq, QS_para);
        Main("dataset/Data_10K.txt", "10K", log, i, QS_seq, QS_para);
        Main("dataset/Data_100K.txt", "100K", log, i, QS_seq, QS_para);

        // 关闭文件以使得修改可见
        log->~Log();        
    }
}