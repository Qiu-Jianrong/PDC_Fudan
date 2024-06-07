# include<bits/stdc++.h>
# include<chrono>

using namespace std;
using namespace chrono;


// 读入数据，待排序的数组
void readData(const string& filename, vector<int>& nums){
    nums.clear();
    string data;

    ifstream file;
    file.open(filename, ios::in);
    while(file){
        file >> data;
        nums.push_back(stoi(data));
    }
    file.close();
}

// 读入矩阵数据，本次实验只使用方阵
void readMatrix(const string& filename, vector<vector<float>>& nums, int n){
    ifstream infile;
    infile.open(filename, ios::in);

    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            float num;
            infile >> num;
            nums[i][j] = num;
        }
    }
}

// 正确性检验，差分测试
bool CorrectnessTest(vector<int>& nums1, vector<int>& nums2){
    int n = nums1.size();
    for(int i = 0; i < n; ++i){
        if (nums1[i] != nums2[i])
            return false;
    }
    return true;
}

// 矩阵相乘正确性检验，差分测试
bool MatrixCorrectness(const vector<vector<float>>& C1, const vector<vector<float>>& C2, int n){

    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            if(fabs(C1[i][j] - C2[i][j]) > 1e-6)
                return false;
        }
    }
    return true;
}
