# include<bits/stdc++.h>
using namespace std;

int main(){
    srand(time(NULL));
    ofstream file;
    file.open("Data_1024.txt", ios::out);
    for(int i = 0; i < 32*32; ++i){
        for(int j = 0; j < 32*32; ++j){
            // 0~2之间的浮点数矩阵，三位有效数字
            file << setprecision(3) << ((float)(rand() % 1024)) / 512 << ' '; 
        }
        file << endl;
    }
    file.close();
}