#pragma once
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <iterator>
#include <cassert>
#include <string>
#include <cmath>
#include <functional>
#include <memory>
#include <random>
#include "Init.hpp"

class MRNA {
public:
    MRNA(const std::string& name, int time,bool cycle) : name(name), time(time), cellCycle(cycle) {}
    std::string getName() const { return name; }
    int getTime() const { return time; }
    bool isCycleMRNA() const { return cellCycle; }
private:
    std::string name;
    int time;
    bool cellCycle = false;
};

extern double mrna_half_life;
extern double mrna_amplification;
extern double mrna_degradation;

class Cell {

public:
    Cell(std::string state ) : id(nextId), cloneID(nextId++), state(state),cycleState(0) {}
    Cell(int cloneID, std::string state) : id(nextId++), cloneID(cloneID), state(state), cycleState(0)   {}
    Cell(int cloneID, std::string state,int cycleState, std::vector<MRNA> mrnaList,bool simMRNA,int cycleCounter) : 
    id(nextId++), cloneID(cloneID), state(state), cycleState(cycleState), mrnaList(mrnaList), simMRNA(simMRNA), cycleCounter(cycleCounter)   {}



    static std::mt19937 rng;
    static std::uniform_real_distribution<double> dist;
    static int nextId;
    static int cell_cycle_duration;



    int getId() const { return id; }
    int getCloneId() const { return cloneID; }
    std::string getState() const { return state; }
    int getCycleState() const { return cycleState; } 
    void incrementCycleCounter() { cycleCounter++; }
    int getCycleCounter() const { return cycleCounter; }
    const std::vector<MRNA>& getMRNAList() const { 
        return mrnaList; 
    }   
    void simulateMRNA(bool sim) { 
        simMRNA = sim; 
    }      
    int getMRNAListSize() const { 
        return mrnaList.size(); 
    }       

    bool maybeStartCellCycle(float probCycle) {        
        if (dist(rng) < probCycle) {
            return true;
        }
        return false;
    }

    Cell* clone() {
        // Create a new cell with the same properties
        Cell* newCell = new Cell(cloneID,state);
        newCell->cycleState = cycleState; // Copy the cycle state
        newCell->mrnaList = mrnaList; // Copy the mRNA list
        newCell->simMRNA = simMRNA; // Copy the mRNA simulation flag
        return newCell;
    }   

    void simulate(std::map<std::string, CellType>& cellTypes,
        std::vector<Cell>&results,
        std::map<std::string, int>& cellTypeFlux,
        int simulationTime);

    void simulateMRNA(int simulationTime);
    int mRNA_creation_rule(double rate, std::mt19937& rng);
    void simulateCycleMRNA(int simulationTime);
    void degradeMRNA(int simulationTime);

private:
    int id;
    int cloneID;
    std::string state;
    int duplicateCount = 0; // Count of duplicates for this cell 
    int cycleState;
    std::vector<MRNA> mrnaList;
    bool simMRNA = false; // Flag to indicate if mRNA simulation is enabled
    int cycleCounter = 0;

};