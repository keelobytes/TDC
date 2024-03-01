Insertion sort optimized using OpenMP pragma tags 

````````

// Optimized

#include <omp.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>


void parallelInsertionSort(std::vector<int>& arr) {
#pragma omp parallel default(none) shared(arr)
    {
        int n = arr.size(); // Define 'n' within the parallel region

#pragma omp for
        for (int i = 1; i < n; ++i) {
            int key;
            key = arr[i];
            int j = i - 1;

            while (j >= 0 && arr[j] > key) {
                arr[j + 1] = arr[j];
                j = j - 1;
            }

            arr[j + 1] = key;
        }
    }
}


int main() {
    omp_set_num_threads(32);
    const int ARRAY_SIZE = 150000;
    std::vector<int> data;
    data.reserve(ARRAY_SIZE);

    // Parallelize the initialization loop using OpenMP
#pragma omp parallel for default(none) shared(data, ARRAY_SIZE)
    for (int i = 0; i < ARRAY_SIZE; ++i) {
        int random_value;
#pragma omp critical
        {
            random_value = std::rand() % 1000;
        }
#pragma omp critical
        {
            data.push_back(random_value);
        }
    }

    clock_t start = clock();

    parallelInsertionSort(data);

    clock_t end = clock();
    double elapsed_secs = static_cast<double>(end - start) / CLOCKS_PER_SEC;

    std::cout << "Time taken: " << elapsed_secs << " seconds" << std::endl;

    return 0;
}


```
