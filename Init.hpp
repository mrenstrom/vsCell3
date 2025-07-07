#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
struct CellType {
    std::string state;
    std::vector<std::pair<std::string, double>> tranList; // Differentiate type and probability
    float probCycle;
    double probDifferentiate;
};

std::map<std::string, CellType> readCellTypesFromCSV(const std::string& filename);
std::map<std::string, float> readMRNAFromCSV(const std::string& filename);
std::map<std::string, std::map<std::string, float>> buildRnaList();


