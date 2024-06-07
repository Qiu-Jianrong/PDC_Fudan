# include<bits/stdc++.h>
using namespace std;
class Log
{
private:
    fstream file;
public:
    Log(const string& filename);
    ~Log();

    void write(const string& batch, double seq, double para){
        file << "Batch size: " << batch << 
            "\nSequential time: " << seq << " ms\n" <<
            "Parallel time: " << para << " ms\n" << 
            "Speedup: " << seq / para << "\n\n";
    }

    // 写入矩阵乘法结果
    void writeMatrix(int& batch, double seq, double cannon){
        file << "Matrix size: " << batch << " * " << batch <<
        "\nSequential time: " << seq << "ms\n" <<
        // "Cannon Algorithm time: " << cannon << "ms\n" <<
        "Simpl Block time: " << cannon << "ms\n" <<
        "Speedup: " << seq / cannon << "\n\n";
    }

};

Log::Log(const string& filename)
{
    // 只写不读，写前清空
    this->file.open(filename, ios::app);
}

Log::~Log()
{
    this->file.close();
}
