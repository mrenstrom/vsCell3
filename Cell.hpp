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
    // Cell(std::string state ) : id(nextId), cloneID(nextId++), state(state),cycleState(0) {}


    // Cell(int cloneID, std::string state) : id(nextId++), cloneID(cloneID), state(state), cycleState(0)   {}


    Cell(int cloneID, std::string state,int cycleState, std::vector<MRNA> mrnaList,bool simMRNA,bool rnaEnabledForType,int cycleCounter,int duplicateCount = 1) 
         : id(nextId++), 
         cloneID(cloneID), 
         state(state), 
         cycleState(cycleState), 
         mrnaList(mrnaList), 
         simMRNA(simMRNA), 
         rnaEnabledForType(rnaEnabledForType), 
         cycleCounter(cycleCounter), 
         duplicateCount(duplicateCount) {}

    Cell(int cloneID, std::string state,Init& g_init) {
        id = nextId++;
        this->cloneID = cloneID;
        this->state = state;
        cycleState = 0; // Default cycle state
        rnaEnabledForType = g_init.rnaEnabledForType(state); // Check if RNA is enabled for this cell type
        simMRNA = false;
        cycleCounter = 0; // Initialize cycle counter
        duplicateCount = g_init.duplicateCountForType(state); // Initialize duplicate count
        // Initialize mRNA list based on the cell type
    }


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


    void setState(const std::string& newState, Init& g_init) { 
        state = newState;
        if (newState.find("Mature") != std::string::npos) {
            // If the new state is "Mature", disable mRNA simulation
            simMRNA = false; // Disable mRNA simulation for mature cells
            mrnaList.clear(); // Clear mRNA list for mature cells
            rnaEnabledForType = false; // RNA is not enabled for mature cells
        } else {
            bool doRNA = g_init.rnaEnabledForType(newState); // Update RNA enabled flag based on new state
            if (!doRNA) {
                rnaEnabledForType = false; // RNA is not enabled for the new state
                simMRNA = false; // Disable mRNA simulation if RNA is not enabled for the new state
                mrnaList.clear(); // Clear mRNA list if RNA is not enabled
            }
        }   
    }

    const std::vector<MRNA>& getMRNAList() const { 
        return mrnaList; 
    }   
    void setSimulateMRNA(bool sim) { 
        if (sim) {
            if (rnaEnabledForType) {
                // If RNA is enabled for this cell type, allow mRNA simulation
                simMRNA = true; 
            } else {
                // If RNA is not enabled for this cell type, disable mRNA simulation
                simMRNA = false; 
            }
        } else {
            simMRNA = false; // Disable mRNA simulation
        }
    } 

    bool getSimulateMRNA() const { 
        return simMRNA; 
    }   

    bool getRNAEnabledForType() const { 
        return rnaEnabledForType; 
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

    int getDuplicateCount() const {
        return duplicateCount;
    }   

    void addToDuplicateCount(int count) {
        duplicateCount += count;
    }

    void setDuplicateCount(int count) {
        duplicateCount = count;
    }   



    Cell* clone(Init& g_init) const {
        // Create a new cell with the same properties
        Cell* newCell = new Cell(cloneID,state,g_init);
        newCell->cycleState = cycleState; // Copy the cycle state
        newCell->mrnaList = mrnaList; // Copy the mRNA list
        newCell->rnaEnabledForType = rnaEnabledForType; // Copy the RNA enabled flag
        newCell->simMRNA = simMRNA; // Copy the mRNA simulation flag
        newCell->cycleCounter = cycleCounter; // Copy the cycle counter
        newCell->duplicateCount = duplicateCount; // Copy the duplicate count
        // Note: The id will be automatically assigned by the static nextId variable
        // so we don't need to set it here.
        return newCell;
    }   

    void simulate(Init  & g_init,
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
    int duplicateCount = 1; // Count of duplicates for this cell 
    int cycleState = 0; // Current cycle state, 0 means not in cycle
    std::vector<MRNA> mrnaList;
    int cycleCounter = 0;  
    bool rnaEnabledForType = false; // Flag to indicate if RNA is enabled for this cell type
    bool simMRNA = false; // Flag to indicate if mRNA simulation is enabled
};