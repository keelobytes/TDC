Insertion Sort Original 

`````````

// Original

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

void insertionSort(std::vector<int>& arr) {
    int n = arr.size();
    for (int i = 1; i < n; ++i) {
        int key = arr[i];
        int j = i - 1;

        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }

        arr[j + 1] = key;
    }
}

int main() {
    // Generate a large unsorted array
    const int ARRAY_SIZE = 150000;
    std::vector<int> data;
    data.reserve(ARRAY_SIZE);

    std::srand(static_cast<unsigned>(std::time(nullptr)));
    for (int i = 0; i < ARRAY_SIZE; ++i) {
        data.push_back(std::rand() % 1000);
    }

    // Uncomment the line below if you want to see the unsorted array
    // for (int num : data) { std::cout << num << " "; }

    clock_t start = clock();

    insertionSort(data);

    clock_t end = clock();
    double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;

    // Uncomment the line below if you want to see the sorted array
    // for (int num : data) { std::cout << num << " "; }

    std::cout << "Time taken: " << elapsed_secs << " seconds" << std::endl;

    return 0;
}

```````
