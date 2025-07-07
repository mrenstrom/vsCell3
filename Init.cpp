#include "Init.hpp"


std::map<std::string, CellType> readCellTypesFromCSV(const std::string& filename) {

    std::map<std::string, CellType> cellTypeMap;
    std::ifstream file(filename);
    std::string line;


    

    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return cellTypeMap;
    }

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string state, type, probCycleStr, probDiffStr;
        std::string tranStateStr, tranProbStr;

        if (!std::getline(ss, state, ',')) continue;
        if (!std::getline(ss, type, ',')) continue;
        if (type == "C") {
            if (!std::getline(ss, probCycleStr, ',')) {
                std::cerr << "Error: Missing probCycle for state '" << state << "' in line: " << line << std::endl;
                continue;
            }
            if (!std::getline(ss, probDiffStr, ',')) {
                std::cerr << "Error: Missing probDifferentiate for state '" << state << "' in line: " << line << std::endl;
                continue;
            }   
            CellType cell;
            cell.state = state;
            cell.probCycle = std::stof(probCycleStr);
            cell.probDifferentiate = std::stod(probDiffStr);
            cellTypeMap[state] = cell;

        } else if (type == "T") {
            if (!std::getline(ss, tranStateStr, ',')) {
                std::cerr << "Error: Missing tranState for state '" << state << "' in line: " << line << std::endl;
                continue;   
            }
            if (!std::getline(ss, tranProbStr, ',')) {
                std::cerr << "Error: Missing tranProb for state '" << state << "' in line: " << line << std::endl;
                continue;   
            }   
            double tranProb = std::stod(tranProbStr);
            if (cellTypeMap.find(state) == cellTypeMap.end()) {
                std::cerr << "Error: State '" << state << "' not found in cell type map." << std::endl;
                std::cerr << "Line: " << line << std::endl;
                continue;
            }   
            CellType cell = cellTypeMap[state]; // Get existing cell type
            if (cell.state.empty()) {
                std::cerr << "Error: State '" << state << "' not found in cell type map." << std::endl;
                std::cerr << "Line: " << line << std::endl;
                continue;
            }   
            cell.tranList.push_back({tranStateStr, tranProb});
            cellTypeMap[state] = cell; // Update the cell type in the map
        } else {
            std::cerr << "Error: Invalid cell type '" << type << "' in line: " << line << std::endl;
            continue;
        }
    }

    return cellTypeMap;
}
//
//
// read mRNA definitions from CSV file
//
std::map<std::string, float> readMRNAFromCSV(const std::string& filename) {
    std::map<std::string, float> mRNAList;
    std::ifstream file(filename);
    std::string line;
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return mRNAList;
    } 
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string idStr, expression;
        if (!std::getline(ss, idStr, ',')){
            std::cerr << "Error: Missing ID in line: " << line << std::endl;
            continue;       
        };
        if (!std::getline(ss, expression, ',')){
            std::cerr << "Error: Missing expression in line: " << line << std::endl;
            continue;           
        }

        float fExp = std::stof(expression);
        mRNAList[idStr] = fExp; // Store the mRNA with its ID and expression value
    }
    return mRNAList;
}
//
// build mRNA expression list
//
std::map<std::string, std::map<std::string, float>> buildRnaList() {
    std::vector<std::pair<std::string,std::string>> mrnaFileNames = {
        {"BAS","bm6_bas_expression.csv"},
        {"CD14MONO","bm6_cd14mono_expression.csv"},
        {"ER1","bm6_er1_expression.csv"},
        {"ER2","bm6_er2_expression.csv"},
        {"MY1","bm6_my1_expression.csv"},
        {"MY2","bm6_my2_expression.csv"},
        {"PREB","bm6_preB_expression.csv"},
        {"PROB","bm6_proB_expression.csv"},
        {"HSC","bm6_stem_expression.csv"},
        {"S","bm6_s_expression.csv"},
        {"G2M","bm6_g2m_expression.csv"}
    };

    std::map<std::string, std::map<std::string, float>> mRNAExpressionList;
    for (const auto& pair : mrnaFileNames) {
        const std::string& cellType = pair.first;
        const std::string& fileName = pair.second;
        std::string filePath = "/Users/markenstrom/Documents/BM5BasedModel/R_Model_estimate/" + fileName;
        std::map<std::string, float> mRNAList = readMRNAFromCSV(filePath);
        if (!mRNAList.empty()) {
            mRNAExpressionList[cellType] = mRNAList;
        } else {
            std::cerr << "Warning: No mRNA data found for cell type " << cellType << " in file " << filePath << std::endl;
        }
    }
    return mRNAExpressionList;
}

 