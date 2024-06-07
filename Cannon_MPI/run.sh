#!/bin/bash
cd C:\\Users\\Qiu\\Desktop\\Project\\Cannon_MPI

g++ main.cpp -o main -I D:\\MPI\\Include -L D:\\MPI\\Lib\\x64 -lmsmpi
for i in {4,16,64}
do
./main.exe -n $i
done
# ./main.exe -n 2
# read -n 1

rm *.exe