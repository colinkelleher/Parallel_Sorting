#!/bin/bash
 echo ==========================
 echo "MPI Sorting"
 echo ==========================

# ==================================
# Compiling the sorts
# ==================================
 echo "Compiling the sorts - simpleSort, BucketSort, oddEven, MergeSort"
 echo ==========================
 mpicc simpleSort.c -o simpleSortEXE
 mpicc oddEven.c -o oddEvenEXE
 mpicc BucketSort.c -o BucketSortEXE
 mpicc MergeSort.c -o MergeSortEXE -lm

# ==================================
# Running the sorts
# ==================================
echo "Running MERGE SORT with $1 processor(s) ";
mpirun -np $1 MergeSortEXE

echo ==========================
echo "Running SIMPLE SORT with $1 processor(s) ";
mpirun -np $1 simpleSortEXE

echo ==========================
echo "Running ODD-EVEN SORT with $1 processor(s) ";
mpirun -np $1 oddEvenEXE

echo ==========================
echo "Running BUCKET SORT with $1 processor(s) ";
mpirun -np $1 BucketSortEXE

echo ==========================
echo "Sorting Complete!";
# ==================================
# Removing the compiled files after running
# ==================================
rm simpleSortEXE
rm oddEvenEXE
rm BucketSortEXE
rm MergeSortEXE
