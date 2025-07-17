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
    // Add a count for duplicates of this cell type
    int duplicateCount = 0;
    bool rnaEnabled = false; // Flag to indicate if this cell type has RNA
};

//std::map<std::string, CellType> readCellTypesFromCSV(const std::string& filename);
std::map<std::string, float> readMRNAFromCSV(const std::string& filename);
std::map<std::string, std::map<std::string, float>> buildRnaList();

class Init {
public:
    void initializeCellTypes(std::string cellTypeFile) {
        readCellTypesFromCSV(cellTypeFile);
    }
    std::string cellTypeFile;
    std::string initialCellType;
    int InitialCellNumber = 0;
    std::map<std::string, CellType> cellTypes;
    void readCellTypesFromCSV(const std::string& filename);
    bool rnaEnabledForType(std::string cellType) {
        auto it = cellTypes.find(cellType);
        if (it != cellTypes.end()) {
            return it->second.rnaEnabled;
        }
        std::cerr << "Error: Cell type '" << cellType << "' not found in cell types map." << std::endl;
        return false; // Default to false if not found
    }
    int duplicateCountForType(std::string cellType) {
        auto it = cellTypes.find(cellType);
        if (it != cellTypes.end()) {
            return it->second.duplicateCount;
        }
        std::cerr << "Error: Cell type '" << cellType << "' not found in cell types map." << std::endl;
        return 0; // Default to 0 if not found
    }
};