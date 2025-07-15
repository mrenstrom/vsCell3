#include <iostream>
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
#include <thread>
#include <regex>
#include "Fred.hpp"
#include "Cell.hpp"
#include "Init.hpp"
//
// Read mRNA expression data from CSV files
//
std::map<std::string,std::map<std::string, float>> mRNAExpressionList;
Init g_init;
//
// Class Cell init random number generator and distribution
//
std::mt19937 Cell::rng(std::random_device{}());
std::uniform_real_distribution<double> Cell::dist(0.0, 1.0);
int Cell::nextId = 0;
int Cell::cell_cycle_duration = 2; // default 2 days
// should these be CELL static?
double mrna_half_life = 1.3; // Default half-life for mRNA in days
double mrna_amplification = 0.65; // Default amplification factor for mRNA
double mrna_degradation = 0.5; // Default degradation rate for mRNA
//
// output all cell type and mrna information to a csv file
//
void outputCellTypesAndMRNA(const std::vector<std::vector<Cell>>& cellGroups,
                            const std::string& filename) {
    std::cout << "Outputting cell types and mRNA to " << filename << std::endl;

    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Cannot open output file." << std::endl;
        return;
    }
    // Write header
    outFile << "Cell ID";
    // mrna names
    std::vector<std::string> allMrnaNames;
    auto mrnaNames = mRNAExpressionList.begin()->second;
    for (const auto& pair : mrnaNames) {
        outFile << "," << pair.first;
        allMrnaNames.push_back(pair.first);
    }
    mrnaNames =mRNAExpressionList["S"];
    for (const auto& pair : mrnaNames) {
        outFile << "," << pair.first;
        allMrnaNames.push_back(pair.first);
    }
    mrnaNames =mRNAExpressionList["G2M"];
    for (const auto& pair : mrnaNames) {
        outFile << "," << pair.first;
        allMrnaNames.push_back(pair.first); 
    }


    outFile << "\n";

    for (const auto& cells:cellGroups) {
        for (const auto& cell : cells) {
            outFile << cell.getState() << "_"
                    << cell.getId() << "_"
                    << cell.getCloneId() << "_"
                    << cell.getCycleState();
            //
            // collect mRNA sums
            //
            std::map<std::string, int> mrnaCount;
            for (const auto& mrna : cell.getMRNAList()) {
                // Count the occurrences of each mRNA type
                mrnaCount[mrna.getName()]++;
            }
            // Write mRNA counts for this cell
            for (const auto& mrnaName : allMrnaNames) {
                // Check if this mRNA type exists in the cell's mRNA list
                if (mrnaCount.find(mrnaName) != mrnaCount.end()) {
                    outFile << "," << mrnaCount[mrnaName];
                } else {
                    outFile << ",0"; // If not found, write 0
                }
            }
            outFile << "\n"; // New line for the next cell            
        }
    }
    outFile.close();
}

//
// Simulate a range of cells in parallel
// This function simulates a range of cells and stores the results in a vector.
// Each thread will process a chunk of cells, simulating each cell and storing the results.
// If a cell divides, the new cell is added to the results vector.
// The original cell is also added to the results vector, so the results will contain both the
// original cells and any new cells that were created during the simulation.
//
// simulate a range of cells
//
void simulate_range(std::vector<Cell>& cells,
    std::map<std::string, int>& cellTypeFlux,
    int startDay,
    int stopDay) {
    
    //
    // remember flux of cell types
    //
    for (const auto& pair : g_init.cellTypes) {
        cellTypeFlux[pair.first] = 0; // Initialize flux for each cell type
    }
    //
    // simulate each day
    //
    for (int simulationTime = startDay; simulationTime < stopDay; ++simulationTime) {
        //
        // Initialize a vector to store the results of the simulation
        //
        std::vector<Cell> results; 
        for (auto& cell : cells) {
            std::string originalState = cell.getState();
            //
            // simulate the cell, if the cell divides, it will add a new cell to results
            // also if the cell type changes, a new cell will be added to results
            cell.simulate(g_init.cellTypes, results, cellTypeFlux, simulationTime);
        }
        // 
        // add the results (daughter cells) to the cells vector
        //
        //cells.insert(cells.end(), std::make_move_iterator(results.begin()), std::make_move_iterator(results.end()));
        for (auto it = results.begin(); it != results.end(); ) {
            std::string cellState = it->getState();
            // if not mature, add to cells
            if (cellState.find("Mature") == std::string::npos) {
                cells.push_back(std::move(*it)); // Move the element
            } 
            ++it;
        }   
        results.clear();    
        //
        // remove cells marked for erasure
        // and also find out how many cells may be combined in duplicates
        //
        std::map<std::string, int> duplicateCount;       
        for (auto it = cells.begin(); it != cells.end();) {
            std::string cellState = it->getState();
            if (cellState.find("erase") != std::string::npos)  {
                it = cells.erase(it); // Remove the cell
            } else {            
                ++it; // Move to the next cell
            }
            int stateDuplicateCount = g_init.cellTypes[cellState].duplicateCount;
            if (stateDuplicateCount > 1) {
                // Increment the count for this cell type
                duplicateCount[cellState]++;
            }
        }
        //
        // show duplicate counts
        //
        //for (const auto& pair : duplicateCount) {
        //    std::cout << "Cell type: " << pair.first 
        //              << ", Duplicate count: " << pair.second << std::endl;
        //}
    }
}
//
// simulate the cells for a given number of days
//  
void simulateCells(std::vector<std::vector<Cell>>& cellGroups,
                   int cores,
                   int startDay,
                   int stopDay) {

    std::vector<std::thread> threads;
    // 
    // Create a map to store the flux of cell types for each thread
    //
    std::vector<std::map<std::string, int>> cellTypeFlux(cores);
    //
    // Launch threads to simulate each group of cells
    //
    for (size_t i = 0; i < cellGroups.size(); ++i) {
        threads.emplace_back(simulate_range, std::ref(cellGroups[i]), std::ref(cellTypeFlux[i]), startDay, stopDay);
    }
    //
    // Wait for all threads to finish
    //
    for (auto& t : threads) {
        t.join();
    }
    //
    // sum up flux vectors
    //
    std::map<std::string, int> totalCellTypeFlux;
    for (const auto& flux : cellTypeFlux) {
        for (const auto& pair : flux) {
            totalCellTypeFlux[pair.first] += pair.second;
        }
    }   
    for (const auto& pair : totalCellTypeFlux) {
        std::cout << "Cell type: " << pair.first << ", Flux: " << pair.second << std::endl;
    }       
    //
    // Count the total number of cells across all groups
    //
    int totalCells = 0;
    for (const auto& cellGroup : cellGroups) {
        totalCells += cellGroup.size();
    }
    std::cout << "Completed days " << startDay << " to " << stopDay << ". Total cells: " << totalCells << std::endl;
}

void handleResize(int) {
    // Do nothing, just override default handler
}
#include <csignal>
//
// Main function to run the vsCell simulation
//
int main(int argc, char* argv[]) {
    std::signal(SIGWINCH, handleResize); // Ignore or safely handle window resize
    std::cout << "Starting vsCell simulation" << " argc = " << argc << std::endl;
    std::string init_cell_types = "cell_types.csv";


    if (argc > 1) {
        init_cell_types = argv[1];
        std::cout << "Using cell types file: " << init_cell_types << std::endl; 
        //
        // look for other args
        //
        for (int i = 2; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--half-life" && i + 1 < argc) {
                mrna_half_life = std::stod(argv[++i]);
                std::cout << "Using mRNA half-life: " << mrna_half_life << " days" << std::endl;
            } else if (arg == "--amplification" && i + 1 < argc) {
                mrna_amplification = std::stod(argv[++i]);
                std::cout << "Using mRNA amplification factor: " << mrna_amplification << std::endl;
            } else if (arg == "--cell_cycle_duration" && i + 1 < argc) {
                Cell::cell_cycle_duration = std::stoi(argv[++i]);
                std::cout << "Using cell cycle duration: " << Cell::cell_cycle_duration << " days" << std::endl;
            } else {
                std::cerr << "Unknown argument: " << arg << std::endl;
            }
        }
    } else {
        std::cout << "Using default cell types file: " << init_cell_types << std::endl; 
    }   
    // 
    // determine the number of available hardware threads
    //
    unsigned int cores = std::thread::hardware_concurrency();
    std::cout << "Available hardware threads: " << cores << std::endl;

    //cores = 1; // Force single-threaded execution for debugging

    //
    // Read cell types from CSV file
    //
    g_init.initializeCellTypes(init_cell_types);
    std::cout << "Number of cell types read: " << g_init.cellTypes.size() << std::endl;
    //
    // Read mRNA expression data from CSV files
    //
    mRNAExpressionList = buildRnaList();
    std::cout << "Number of mRNA types read: " << mRNAExpressionList.size() << std::endl;   
    //
    // init cells into mutliple vectors for parallel processing
    //    
    std::vector<std::vector<Cell>> cellGroups(cores);
    int chunk_size = g_init.InitialCellNumber / cores;
    //
    // starting cell type is first line in cellTypes

    if (g_init.cellTypes.empty()) {
        std::cerr << "Error: No cell types found in the CSV file." << std::endl;
        return 1; // Exit with an error code
    } 
    std::cout << "Initial cell type: " << g_init.initialCellType<< std::endl;
    //
    // Create initial cells and distribute them across the available cores
    //
    for (int i = 0; i < cores; ++i) {
        size_t start = i * chunk_size;
        size_t end = (i == cores - 1) ? g_init.InitialCellNumber : start + chunk_size;
        //std::cout << "Creating cells from " << start << " to " << end << std::endl;
        for (size_t j = start; j < end; ++j) {
            // Create a new cell with the initial state and add it to the group
            cellGroups[i].emplace_back(j, g_init.initialCellType);
        }
    }
    std::cout << "Initial number of cells: " << g_init.InitialCellNumber << std::endl;
    //
    // prepare to run the simulation
    //
    std::string command;
    int day = 0;
    int param1 = 0;
    std::string script = "";
    std::ifstream scriptFile;
    while (true)
    {
        //
        // if running a script, read the script file otherwise read user input  
        //
        std::string line;
        if (!script.empty()) {
            if (scriptFile) {
                if (!std::getline(scriptFile, line)) {
                    std::cout << "Script complete." << std::endl;
                    script = "";
                    continue;
                }
            }
        } else {
            std::cout << "Enter command: ";
            std::getline(std::cin, line);  // Read the full line

            std::cout << "Input length = " << line.length() << std::endl;
            if (line.empty()) {
                std::cout << "No input provided. Please try again." << std::endl;
                continue; // Skip to the next iteration
            }
        }

        std::istringstream iss(line);
        std::string command;
        std::string arg;

        iss >> command >> arg;
        std::cout << "Command: " << command << ", Argument: " << arg << std::endl;
        if (command.empty()) {
            //std::cout << "No command entered. Please try again." << std::endl;
            continue; // Skip to the next iteration
        } else if ((command == "q") || (command == "quit")) {
            std::cout << "Exiting." << std::endl;
            return 0; // Exit the program
        } else if (command == "z") {
            //
            // run a command script
            // 
            if (arg.empty()) {
                script = "script.txt"; // Default script file
            } else {
                script = arg; // Use the provided script file
            }
            scriptFile.open(script);
            if (!scriptFile) {
                std::cerr << "Error opening script file: " << script << std::endl;
                script = ""; // Reset script if file cannot be opened
                continue; // Skip to the next iteration
            }
        } else if ((command == "s") || (command == "sim")) {
            //
            // simulate
            //
            if (arg.empty()) {
                std::cout << "Please provide number of days to run." << std::endl;
                continue; // Skip to the next iteration
            }
            param1 = std::stoi(arg); // Convert the argument to an integer
            simulateCells(cellGroups, cores, day, day + param1);
            day += param1; // Update the current day
            std::cout << "Simulation complete for " << param1 << " days." << std::endl;
        } else if (command == "o") {
            std::stringstream ss;
            ss << "/Users/markenstrom/Documents/BM5BasedModel/c/vsCell3/out/simOutput_" << mrna_amplification << "_" << mrna_half_life << ".csv";
            std::string filename = ss.str();
            //
            if (!arg.empty()) {
                filename = arg; // Use the provided filename if given
            }
            //
            // Output cell types and mRNA counts to CSV file
            //
            outputCellTypesAndMRNA(cellGroups, filename);
            continue; // Skip to the next iteration
        } else if (command == "r") {
            //
            // turn on mRNA simulation
            //
            for (auto& cellV : cellGroups) {
                for (auto& cell : cellV) {
                    cell.simulateMRNA(true); // Enable mRNA simulation
                }
            }
            std::cout << "mRNA simulation enabled." << std::endl;
        } else if (command == "t") {
            // show cell types
            std::map<std::string, int> cellTypeCount;
            for (const auto& cellV : cellGroups) {
                for (const auto& cell : cellV) {
                    cellTypeCount[cell.getState()]++;
                }
            }            
            std::cout << "Cell types:" << std::endl;
            for (const auto& pair : cellTypeCount) {
                std::cout << pair.first << ": " << pair.second << std::endl;
                // for this cell type, how many are in cell cycle
                int inCycle = 0;
                for (const auto& cellV : cellGroups) {
                    for (const auto& cell : cellV) {
                        if (cell.getState() == pair.first && cell.getCycleState() > 0) {
                            inCycle++;
                        }
                    }
                }
                std::cout << "  In Cycle: " << inCycle << " r = " << inCycle / static_cast<double>(pair.second) << std::endl;
            }
        } else if (command == "e") {
            if (arg.empty()) {
                // list cellIDs
                for (const auto& cellV : cellGroups) {
                    for (const auto& cell : cellV) {
                        std::cout << "Cell ID: " << cell.getId() 
                                  << ", Clone ID: " << cell.getCloneId() 
                                  << ", State: " << cell.getState() 
                                  << ", Cycle State: " << cell.getCycleState() 
                                  << ", mRNA Count: " << cell.getMRNAListSize()
                                  << ", Cycle Count: " << cell.getCycleCounter()
                                  << std::endl;
                    }
                }
                continue; // Skip to the next iteration
            }
            // check if arg is in CellTypes
            if (g_init.cellTypes.find(arg) != g_init.cellTypes.end()) {
                // list HSC cells
                for (const auto& cellV : cellGroups) {
                    for (const auto& cell : cellV) {
                        if (cell.getState() == arg) {
                            std::cout << "Cell ID: " << cell.getId() 
                                      << ", Clone ID: " << cell.getCloneId() 
                                      << ", State: " << cell.getState() 
                                      << ", Cycle State: " << cell.getCycleState() 
                                      << ", mRNA Count: " << cell.getMRNAListSize() 
                                      << ", Cycle Count: " << cell.getCycleCounter()
                                      << std::endl;
                        }
                    }
                }
                continue; // Skip to the next iteration
            }
            // arg is a cellID
            param1 = std::stoi(arg); // Convert the argument to an integer
            if (param1 <= 0) {
                std::cout << "Invalid cellID. Please provide a positive integer." << std::endl;
                continue; // Skip to the next iteration
            } 
            std::cout << "Outputting cell information for Cell ID: " << param1 << std::endl;
            for (const auto& cellV : cellGroups) {
                for (const auto& cell : cellV) {
                    if (cell.getId() == param1) {
                        std::cout << "Cell ID: " << cell.getId() 
                                  << ", Clone ID: " << cell.getCloneId() 
                                  << ", State: " << cell.getState() 
                                  << ", Cycle State: " << cell.getCycleState() 
                                  << ", mRNA Count: " << cell.getMRNAListSize() 
                                  << std::endl;
                        const auto& mrnaList = cell.getMRNAList();
                        if (!mrnaList.empty()) {
                            std::cout << "mRNA List:" << std::endl;
                            for (const auto& mrna : mrnaList) {
                                std::cout << "  - " << mrna.getName() 
                                          << " (Time: " << mrna.getTime() << ")" 
                                          << std::endl;
                            }
                        } else {
                            std::cout << "No mRNA present." << std::endl;
                        }
                        break; // Exit the loop after finding the cell
                    }   
                }
            } 
        } else if (command == "b") {
            // list cellIDs
            std::vector<Cell> tmp;
            int totalCells = 0;
            for (auto& cellV : cellGroups) {
                std::cout   << "Cell Group size : " << cellV.size() << std::endl;
                totalCells += cellV.size();
                for (auto& cell : cellV) {
                    tmp.push_back(cell);
                }
                cellV.clear(); // Clear the cell group after moving cells
            }
            std::cout << "Total Cells: " << totalCells << std::endl;
            int chunk_size = totalCells / cores;
            for (int i = 0; i < cores; ++i) {
                size_t start = i * chunk_size;
                size_t end = (i == cores - 1) ? totalCells : start + chunk_size;
                // move cells from tmp to the cellGroups
                for (size_t j = start; j < end; ++j)
                {
                    cellGroups[i].push_back(tmp[j]);
                }
            }
            tmp.clear(); // Clear the temporary vector after moving cells

            //
            // check
            //
            for (auto& cellV : cellGroups) {
                std::cout   << "Cell Group size : " << cellV.size() << std::endl;
            }
        
            continue; // Skip to the next iteration        
        } else if (command == "c") {
            // clone diversity
            std::map<int, int> cloneCount;
            for (const auto& cellV : cellGroups) {
                for (const auto& cell : cellV) {
                    cloneCount[cell.getCloneId()]++;
                }
            }
            std::cout << "Clone diversity:" << std::endl;
            // sort clones by number of cells
            std::vector<std::pair<int, int>> sortedCloneCount(cloneCount.begin(), cloneCount.end());
            std::sort(sortedCloneCount.begin(), sortedCloneCount.end(),
                      [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                          return a.second > b.second; // Sort by count in descending order
                      });
            // output sorted clone counts
            for (const auto& pair : sortedCloneCount) {
                std::cout << "Clone ID: " << pair.first << ", Count: " << pair.second << std::endl;
            }
            std::cout << "Total clones: " << cloneCount.size() << std::endl;
        } else if (command == "i") {
            // reinitialize cells
            int numCells = 1000;
            if (!arg.empty()) {
                numCells = std::stoi(arg); // Convert the argument to an integer
            }
            for (auto& cellV : cellGroups) {
                cellV.clear(); // Clear the cell group before reinitializing
            }
            std::cout << "Reinitializing cells with " << numCells << " cells." << std::endl;
            // Reinitialize cells into multiple vectors for parallel processing   
            int chunk_size = numCells / cores;
            for (int i = 0; i < cores; ++i) {
                size_t start = i * chunk_size;
                size_t end = (i == cores - 1) ? numCells : start + chunk_size;
                std::cout << "Creating cells from " << start << " to " << end << std::endl;
                for (size_t j = start; j < end; ++j) {
                    // Create a new cell with the initial state "HSC" and add it to the group
                    cellGroups[i].emplace_back(j, "HSC");
                }
            }
            std::cout << "Initial number of cells: " << numCells << std::endl;
            day = 0;
        } else if (command == "m") {
            // get mrna count
            for (auto& cellV : cellGroups) {
                for (auto& cell : cellV) {
                    std::cout << "Cell ID: " << cell.getId() << ", mRNA Count: " << cell.getMRNAListSize() << std::endl;
                }
            }
        } else if (command == "h") {
            // half_life
            if (arg.empty()) {
                std::cout << "Current mRNA half-life: " << mrna_half_life << " days." << std::endl;
            } else {
                try {
                    mrna_half_life = std::stod(arg); // Convert the argument to a double
                    std::cout << "Updated mRNA half-life to: " << mrna_half_life << " days." << std::endl;
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Invalid argument for half-life. Please provide a valid number." << std::endl;
                }
            }
        } else if (command == "a") {
            // amplification
            if (arg.empty()) {
                std::cout << "Current mRNA amplification factor: " << mrna_amplification << "." << std::endl;
            } else {
                try {
                    mrna_amplification = std::stod(arg); // Convert the argument to a double
                    std::cout << "Updated mRNA amplification factor to: " << mrna_amplification << "." << std::endl;
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Invalid argument for amplification factor. Please provide a valid number." << std::endl;
                }
            }
        } else if (command == "l") {
            // cell cycle duration
            if (arg.empty()) {
                std::cout << "Current cell cycle duration: " << Cell::cell_cycle_duration << " days." << std::endl;
            } else {
                try {
                    Cell::cell_cycle_duration = std::stoi(arg); // Convert the argument to an integer
                    std::cout << "Updated cell cycle duration to: " << Cell::cell_cycle_duration << " days." << std::endl;
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Invalid argument for cell cycle duration. Please provide a valid integer." << std::endl;
                } 
            }
        } else {
            std::cout << "Unknown command. Please enter 's' to start simulation, 'o' to output cell types and mRNA, 'e' to exit, or 'd' to degrade mRNA." << std::endl;
        }
        // Reset the argument for the next command
        arg.clear();
    }




    // If we reach here, it means the user has chosen to exit

    return 0;
}