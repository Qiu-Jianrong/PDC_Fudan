# include<bits/stdc++.h>
using namespace std;

int main(){
    srand(time(NULL));
    ofstream file;
    file.open("Data_100K.txt", ios::out);
    for(int i = 0; i < 1024 * 100; ++i){
        file << rand() % INT_MAX << " ";
    }
    file.close();
}