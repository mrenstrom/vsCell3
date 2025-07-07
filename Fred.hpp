#include <string>
#include <iostream>
#include <thread>   
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>
#include <map>

class Fred {
public:
    Fred() {
        // Constructor code here
        std::cout << "Fred object created." << std::endl;
    }
    ~Fred() {
        // Destructor code here
        std::cout << "Fred object destroyed." << std::endl;
    }
    void doSomething();

};