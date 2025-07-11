#include "Cell.hpp"
extern std::map<std::string,std::map<std::string, float>> mRNAExpressionList;
//
// returns true if daughter cell was created
//
bool Cell::simulate(std::map<std::string, CellType>& cellTypes, std::vector<Cell>& results, int simulationTime) {
    CellType cellDef = cellTypes[state];
    if (cellDef.state.empty()) {
        std::cerr << "Error: Cell type '" << state << "' not found in cell types map." << std::endl;
        return false; // Return false if the cell type is not found
    }       
    //std::cout << "Cell " << id << " state " << state << "

    bool daughterCreated = false;
    degradeMRNA(simulationTime); // Degrade mRNA molecules at the start of the simulation step

    if (cycleState == 0) {
        // Check if the cell can start a cell cycle
        if (maybeStartCellCycle(cellDef.probCycle)) {  
            cycleState = 1; // Start cell cycle
        }
    } else {
        cycleState++;
        if (cycleState >= Cell::cell_cycle_duration) {
            incrementCycleCounter();
            // Reset the cell state
            cycleState = 0;
            daughterCreated = true; // Indicate that a daughter cell was created
            results.emplace_back(cloneID, state); // Create a new daughter cell with the same cloneID and state
            auto daughterCell = &results.back(); // Reference to the newly added Cell
            daughterCell->cycleState = 0;
            daughterCell->mrnaList = mrnaList; // Copy the mRNA list
            daughterCell->simMRNA = simMRNA; // Copy the mRNA simulation flag
            daughterCell->cycleCounter = cycleCounter; // Copy the cycle counter
            //
            // check both mother cell and daughter cell for differentiation
            //
            if (dist(rng) < cellDef.probDifferentiate) {
                // Randomly choose a differentiation type
                double randValue = dist(rng);
                double cumulativeProb = 0.0;
                for (const auto& tran : cellDef.tranList) {
                    cumulativeProb += tran.second;
                    if (randValue < cumulativeProb) {
                        daughterCell->state = tran.first; // Change state to the differentiated type
                        //std::cout << "daughter Cell " << id << " differentiated to " << daughterCell->state << std::endl;
                        break;
                    }
                }
            } 
            //
            // mother cell
            //
            if (dist(rng) < cellDef.probDifferentiate) {
                // Randomly choose a differentiation type
                double randValue = dist(rng);
                double cumulativeProb = 0.0;
                for (const auto& tran : cellDef.tranList) {
                    cumulativeProb += tran.second;
                    if (randValue < cumulativeProb) {
                        state = tran.first; // Change state to the differentiated type
                        //std::cout << "mother Cell " << id << " differentiated to " << state << std::endl;
                        break;
                    }
                }
            }

        }
    }
    simulateMRNA(simulationTime); // Simulate mRNA if enabled
    simulateCycleMRNA(simulationTime); // Simulate mRNA if in cell cycle
    if (daughterCreated) {
        // If the daughter cell was created during the cycle, simulate mRNA for it
        auto daughterCell = &results.back(); // Reference to the newly added Cell
        daughterCell->simulateMRNA(simulationTime);
        daughterCell->simulateCycleMRNA(simulationTime);
    }
    //
    // check mrnaList for the mother cell
    //
    if ((simMRNA) && (mrnaList.empty()) && (state != "Mature")) {
        // If mRNA simulation is enabled and the mother cell has no mRNA molecules, print  
        std::cerr << "Error: No mRNA molecules found for cell " << id << " state " << state << std::endl;
    }

    //
    // add daughter cell to mRNA cell list
    //
    return daughterCreated; // Return true if a daughter cell was created
}
//
// Simulate mRNA for the cell
//

int Cell::mRNA_creation_rule(double rate, std::mt19937& rng) {
    std::poisson_distribution<int> dist(rate);
    return floor(dist(rng));
}

void Cell::simulateMRNA(int simulationTime) {
    int created = 0; // Counter for created mRNA molecules
    if (!simMRNA) return; // If mRNA simulation is not enabled, do nothing
    if (state == "Mature") return; // If the cell is mature, do not simulate mRNA
    // Simulate mRNA creation based on the cell's state
    auto expList = mRNAExpressionList[state];
    // make sure state is in the mRNAExpressionList
    if (expList.empty()) {
        std::cerr << "Error: No mRNA expression data for state '" << state << "'." << std::endl;
        return;
    }
    double mrnaScale = mrna_amplification; // Scale factor for mRNA expression, adjust as needed
    for (const auto& rna : expList) {
        std::string mrnaName = rna.first;
        float mrnaExpression = rna.second * mrnaScale; // adjust to days
        //
        // is a new mRNA produced?
        //
        int molecules = mRNA_creation_rule(mrnaExpression, rng);
        created += molecules; // Increment the created counter
        if (molecules <= 0) continue; // If no mRNA is produced,
        // skip to the next one
        //        // Create the mRNA molecules and add them to the list
        for (int i = 0; i < molecules; ++i) {
            // Create a new mRNA object with the name and simulation time
            mrnaList.emplace_back(mrnaName, simulationTime,false); // Add the new mRNA to the list
        }
    }
    //std::cout << "Cell " << id << " created " << created << " mRNA molecules at time " << simulationTime << std::endl;
}   
//
// cell cycle mRNA turns on fast
//
void Cell::simulateCycleMRNA(int simulationTime) {
    if (!simMRNA) return; // If mRNA simulation is not enabled, do nothing
    if (state == "Mature") return; // If the cell is mature, do not simulate mRNA
    if (cycleState == 0) return;
    std::map<std::string, float> expList;
    //std::cout << "Cell " << id << " cycleState " << cycleState << std::endl;
    if (cycleState < Cell::cell_cycle_duration / 2) {   
        // S phase
        expList = mRNAExpressionList["S"];
    } else {
        expList = mRNAExpressionList["G2M"];
    } 
    if (expList.empty()) {
        std::cerr << "Error: No mRNA expression data for state cell cycle state." << std::endl;
        return;
    }

    double mrnaScale = mrna_amplification; // Scale factor for mRNA expression, adjust as needed
    for (const auto& rna : expList) {
        std::string mrnaName = rna.first;
        float mrnaExpression = rna.second * mrnaScale * 1.0; // faster creation
        //
        // is a new mRNA produced?
        //
        int molecules = mRNA_creation_rule(mrnaExpression, rng);
        if (molecules <= 0) continue; // If no mRNA is produced,
        // skip to the next one
        //        // Create the mRNA molecules and add them to the list
        for (int i = 0; i < molecules; ++i) {
            // Create a new mRNA object with the name and simulation time
            mrnaList.emplace_back(mrnaName, simulationTime,true); // Add the new mRNA to the list
        }
    }
}
//
//
// Degrade mRNA molecules based on their age and half-life
//
// Returns true if the mRNA molecule survives, false if it degrades
//
bool mRNA_survives(double age_days, double half_life_days, std::mt19937& rng) {
    double k = std::log(2.0) / half_life_days;
    double survival_prob = std::exp(-k * age_days);

    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double p1 = dist(rng);
    return p1 < survival_prob;
}



void Cell::degradeMRNA(int simulationTime) {
    if (!simMRNA) return; // If mRNA simulation is not enabled, do nothing
    // Degrade mRNA molecules based on their expression rate
    double half_life_days = mrna_half_life; // Set the half-life of mRNA molecules in days
    for (auto it = mrnaList.begin(); it != mrnaList.end();) {
        // if (it->isCycleMRNA()) {
        //    if (simulationTime - it->getTime() > 1) {
        //         it = mrnaList.erase(it); // Remove the mRNA molecule if its time is less than the current simulation time
        //    } else {
        //        ++it; // Move to the next mRNA molecule
        //    }           
        // } else 
        if (mRNA_survives(simulationTime - it->getTime(), half_life_days, rng)) {
            it = mrnaList.erase(it); // Remove the mRNA molecule if its time is less than the current simulation time
        } else {
            ++it; // Move to the next mRNA molecule
        }
    }
}   