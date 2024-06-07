#!/bin/bash
cd C:\\Users\\Qiu\\Desktop\\Project\\PSRS_MPI

g++ main.cpp -o main -I D:\\MPI\\Include -L D:\\MPI\\Lib\\x64 -lmsmpi
for ((i = 2; i < 17; i += 2))
do
./main.exe -n $i
done
# ./main.exe -n 2
# read -n 1

rm *.exe