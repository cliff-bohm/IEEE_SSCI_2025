//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#include "WireCubeBrain.h"

std::shared_ptr<ParameterLink<int>> WireCubeBrain::xPL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-size_x", 10, "width of brain matrix");
std::shared_ptr<ParameterLink<int>> WireCubeBrain::yPL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-size_y", 10, "height of brain matrix");
std::shared_ptr<ParameterLink<int>> WireCubeBrain::zPL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-size_z", 10, "depth of brain matrix");

std::shared_ptr<ParameterLink<int>> WireCubeBrain::decayDurationPL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-decayDuration", 1, "how long decay takes");

std::shared_ptr<ParameterLink<int>> WireCubeBrain::stepsPerUpdatePL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-stepsPerUpdate", 1, "number of matrix steps run for each brain update");

std::shared_ptr<ParameterLink<bool>> WireCubeBrain::clearBetweenUpdatesPL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-clearBetweenUpdates", true, "is the brain cleared at the start of each update?");

std::shared_ptr<ParameterLink<double>> WireCubeBrain::initialFillRatioPL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-initialFillRatio", .1, "ratio of cells set to wire in initialized brains");

std::shared_ptr<ParameterLink<int>> WireCubeBrain::nrRecurrentValuesPL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-nrRecurrentValues", 8, "number of recurrent values");

std::shared_ptr<ParameterLink<double>> WireCubeBrain::mutationRate_addWirePL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-mutationRate_addWire", .001, "rate per empty cell to change that cell to wire");
std::shared_ptr<ParameterLink<double>> WireCubeBrain::mutationRate_eraseWirePL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-mutationRate_eraseWire", .01, "rate per wire cell to change that cell to empty");
std::shared_ptr<ParameterLink<double>> WireCubeBrain::mutationRate_alterInputPL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-mutationRate_rewireInput", .01, "rate per input to change which cell input connects to");
std::shared_ptr<ParameterLink<double>> WireCubeBrain::mutationRate_alterOutputPL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-mutationRate_rewireOutput", .01, "rate per output to change which cell output corrects to");

std::shared_ptr<ParameterLink<bool>> WireCubeBrain::optimizeBrainPL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-optimizeBrain", true, "if true, brain will be optimized (e.g., skip unconnected wire) before updates are run\n"
    "optimization should have no effect aside from speed of run");

std::shared_ptr<ParameterLink<int>> WireCubeBrain::saveBrainBehaviorPL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-saveBrainBehavior", 0, "if true, brain behavior will be saved for visulization. This may generate a large file!");

std::shared_ptr<ParameterLink<int>> WireCubeBrain::outputModePL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-outputMode", 0, "if 0, brain outputs are generated from final state of T1_cells\n"
    "if 1, brain outputs are average of T1_cells");

std::shared_ptr<ParameterLink<int>> WireCubeBrain::chargeMinPL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-chargeMin", 1, "a wire must have at least this number of charged neghbiors to become charge");
std::shared_ptr<ParameterLink<int>> WireCubeBrain::chargeMaxPL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-chargeMax", 2, "a wire must have at most this number of charged neghbiors to become charge");

std::shared_ptr<ParameterLink<int>> WireCubeBrain::T0LayersPL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-T0Layers", -1, "which layers (up from y = 0) can have T0 connections(input+hidden), -1 indicates all layers");
std::shared_ptr<ParameterLink<int>> WireCubeBrain::T1LayersPL =
Parameters::register_parameter("BRAIN_WIRE_CUBE-T1Layers", -1, "which layers (down from y = max) can have T1 connections(output+hidden), -1 indicates all layers");


// index to location
std::vector<int> i2l(int index, int xSize, int ySize, int zSize) {
    // x = left/right (width)
    // y = up/down (height)
    // z = in/out (depth)
    std::vector<int> loc(3);

    int xzSize = xSize * zSize;     // size of one layer
    loc[1] = index / xzSize;        // y-coordinate - layer in matrix
    index %= xzSize;                // index in this layer
    loc[2] = index / xSize;         // z-coordinate - row in this layer
    loc[0] = index % xSize;         // x-coordinate - position in this row
    return loc;
}

// location to index (ints)
int l2i(int x, int y, int z, int xSize, int ySize, int zSize) {
    return x + y * (xSize * zSize) + z * xSize;
}

// location to index (vector)
int l2i(std::vector<int> loc, int xSize, int ySize, int zSize) {
    return loc[0] + loc[1] * (xSize * zSize) + loc[2] * xSize;
}

void WireCubeBrain::reportFull() {
    std::cout << "-- Full Report --" << std::endl;
    for (auto i : matrix_wire) {
        std::cout << i;
    }
    auto tempLst = wireCells;
    std::sort(tempLst.begin(), tempLst.end());
    std::cout << std::endl << "wireCells.size() = " << wireCells.size() << std::endl;
    for (auto c : tempLst) {
        std::cout << c << " ";
    }
    std::cout << std::endl << "T0: ";
    for (auto c : T0_cells) {
        std::cout << c << " ";
    }
    std::cout << std::endl << "T1: ";
    for (auto c : T1_cells) {
        std::cout << c << " ";
    }
    std::cout << std::endl;
    std::cout << "links:" << std::endl;
    for (auto c : wireCells) {
        std::cout << "  c(" << c << "): ";
        for (auto l : links[c]) {
            std::cout << l << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "-- End Full Report --" << std::endl;
}

void WireCubeBrain::report() {
    std::cout << "-- Report --" << std::endl;
    auto tempLst = wireCells;
    std::sort(tempLst.begin(), tempLst.end());
    std::cout << "wireCells.size() = " << wireCells.size() << std::endl;
    for (auto c : tempLst) {
        std::cout << c << " ";
    }
    std::cout << std::endl << "T0: ";
    for (auto c : T0_cells) {
        std::cout << c << " ";
    }
    std::cout << std::endl << "T1: ";
    for (auto c : T1_cells) {
        std::cout << c << " ";
    }
    std::cout << std::endl;
    std::cout << "-- End Report --" << std::endl;
}

void WireCubeBrain::showBrainActivityState() {
    if (0) { // show state of cells
        int zz = 0;
        for (auto v : matrix_state) {
            if (matrix_wire[zz] == 0) {
                std::cout << " ";
            }
            else if (v == 0) {
                std::cout << "w";
            }
            else if (v == -1) {
                std::cout << "c";
            }
            else {
                std::cout << "d";
            }
            zz++;
        }
        std::cout << std::endl;
        zz = 0;
        for (auto v : matrix_state_next) {
            if (matrix_wire[zz] == 0) {
                std::cout << " ";
            }
            else if (v == 0) {
                std::cout << "w";
            }
            else if (v == -1) {
                std::cout << "c";
            }
            else {
                std::cout << "d";
            }
            zz++;
        }
        std::cout << std::endl;
    }
}

void WireCubeBrain::makeLinks() {
    if (1) { // 26 way
        // set up links lists (lists for each wire index with wire neighbors, i.e., cells we will need to check)
        links.clear(); // clear out all links lists
        links.resize(cellCount, std::vector<int>(0)); // clear out all links lists
        std::vector<int> offsetLoc(3);
        int offsetIndex;
        
        fullStateBuffer.clear();

        if (saveBrainBehavior){
            fullStateBuffer.resize(cellCount);
        }
        for (auto index : wireCells) { // for each cell that is wire
            if (saveBrainBehavior) {
                fullStateBuffer[index].ID = index; // set this focal ID in state tracker
                //std::cout << "index: " << index << std::endl;
            }
            auto loc = i2l(index, x, y, z); // location in x,y,z of current cell
            for (int xOff : {-1, 0, 1}) { // for every location around index
                for (int yOff : {-1, 0, 1}) {
                    for (int zOff : {-1, 0, 1}) {
                        if (xOff != 0 || yOff != 0 || zOff != 0) { // don't evaluate current cell
                            offsetLoc[0] = loc[0] + xOff; // make an offset location
                            offsetLoc[1] = loc[1] + yOff;
                            offsetLoc[2] = loc[2] + zOff;
                            if (offsetLoc[0] >= 0 && offsetLoc[0] < x && // test that offset location is not off matrix
                                offsetLoc[1] >= 0 && offsetLoc[1] < y &&
                                offsetLoc[2] >= 0 && offsetLoc[2] < z) {
                                offsetIndex = l2i(offsetLoc, x, y, z); // convert offset location to an index
                                if (matrix_wire[offsetIndex] == 1) { // if offset index location is wire
                                    links[index].push_back(offsetIndex); // add offset index to links for current cell
                                    if (saveBrainBehavior) {
                                        //std::cout << "cell: " << index << "(" << loc[0] << "," << loc[1] << "," << loc[2] << ")" << "  has neighbor : " << offsetIndex << "(" << offsetLoc[0] << "," << offsetLoc[1] << "," << offsetLoc[2] << ")" << std::endl;
                                        fullStateBuffer[index].neighborIDs.push_back(offsetIndex); // set this neighbor ID in state tracker
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (0) { // 14 way
    // set up links lists (lists for each wire index with wire neighbors, i.e., cells we will need to check)
        links.clear(); // clear out all links lists
        links.resize(cellCount, std::vector<int>(0)); // clear out all links lists
        std::vector<int> offsetLoc(3);
        int offsetIndex;
        for (auto index : wireCells) { // for each cell that is wire
            auto loc = i2l(index, x, y, z); // location in x,y,z of current cell
            for (int xOff : {-1, 0, 1}) { // for every location around index
                for (int yOff : {-1, 0, 1}) {
                    for (int zOff : {-1, 0, 1}) {
                        if ((xOff != 0 || yOff != 0 || zOff != 0) && (std::abs(xOff) + std::abs(yOff) + std::abs(zOff) != 3)){ // don't evaluate current cell
                            offsetLoc[0] = loc[0] + xOff; // make an offset location
                            offsetLoc[1] = loc[1] + yOff;
                            offsetLoc[2] = loc[2] + zOff;
                            if (offsetLoc[0] >= 0 && offsetLoc[0] < x && // test that offset location is not off matrix
                                offsetLoc[1] >= 0 && offsetLoc[1] < y &&
                                offsetLoc[2] >= 0 && offsetLoc[2] < z) {
                                offsetIndex = l2i(offsetLoc, x, y, z); // convert offset location to an index
                                if (matrix_wire[offsetIndex] == 1) { // if offset index location is wire
                                    links[index].push_back(offsetIndex); // add offset index to links for current cell
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    
    if (0) { // 6 way
             // set up links lists (lists for each wire index with wire neighbors, i.e., cells we will need to check)
        links.clear(); // clear out all links lists
        links.resize(cellCount, std::vector<int>(0)); // clear out all links lists
        std::vector<int> offsetLoc(3);
        int offsetIndex;
        std::vector<std::vector<int>> sixWay;
        sixWay.push_back({ 1, 0, 0 });
        sixWay.push_back({ -1, 0, 0 });
        sixWay.push_back({ 0, 1, 0 });
        sixWay.push_back({ 0, -1, 0 });
        sixWay.push_back({ 0, 0, 1 });
        sixWay.push_back({ 0, 0, -1 });
        for (auto index : wireCells) { // for each cell that is wire
            auto loc = i2l(index, x, y, z); // location in x,y,z of current cell
            for (std::vector<int> offLoc : sixWay) { // for every location around index
                offsetLoc[0] = loc[0] + offLoc[0]; // make an offset location
                offsetLoc[1] = loc[1] + offLoc[1];
                offsetLoc[2] = loc[2] + offLoc[2];
                if (offsetLoc[0] >= 0 && offsetLoc[0] < x && // test that offset location is not off matrix
                    offsetLoc[1] >= 0 && offsetLoc[1] < y &&
                    offsetLoc[2] >= 0 && offsetLoc[2] < z) {
                    offsetIndex = l2i(offsetLoc, x, y, z); // convert offset location to an index
                    if (matrix_wire[offsetIndex] == 1) { // if offset index location is wire
                        links[index].push_back(offsetIndex); // add offset index to links for current cell
                    }
                }
            }
        }
    }

    if (optimizeBrain) { // optimize        
        wireCellsActive.clear();
        for (auto testIndex : T1_cells) { // for each output and T+1 reccurent
            if (std::find(wireCells.begin(), wireCells.end(), testIndex) != wireCells.end()) { // make sure it's wire, i.e., in wireCells
                wireCellsActive.push_back(testIndex); // add this location to list to check, these are linked directly to T+1
            }
        }
        size_t i = 0;
        while (i < wireCellsActive.size()) { // while there are still wires to check
            auto testIndex = wireCellsActive[i];
            for (auto link : links[testIndex]) {
                //std::cout << "cell: " << testIndex << "  has neighbor : " << link << std::endl;
                if (std::find(wireCellsActive.begin(), wireCellsActive.end(), link) == wireCellsActive.end()) {
                    // if not in wireCellsActive
                    wireCellsActive.push_back(link);
                }
            }
            i += 1;
        }
    }
    else {
        wireCellsActive = wireCells;
    }
}

// setup for first call, this is called on every init, but importantly when creating the template brain\
// and when we want an empty brain to copy offspring into
WireCubeBrain::WireCubeBrain(int ins, int outs, std::shared_ptr<ParametersTable> PT)
    : AbstractBrain(ins, outs, PT) {

    nrInputValues = ins;
    nrOutputValues = outs;
    
    x = xPL->get(PT);
    y = yPL->get(PT);
    z = zPL->get(PT);
    cellCount = x * y * z;

    decayDuration = decayDurationPL->get(PT);
    stepsPerUpdate = stepsPerUpdatePL->get(PT);
    clearBetweenUpdates = clearBetweenUpdatesPL->get(PT);
    initialFillRatio = initialFillRatioPL->get(PT);
    nrRecurrentValues = nrRecurrentValuesPL->get(PT);

    T0_size = nrInputValues + nrRecurrentValues;
    T1_size = nrOutputValues + nrRecurrentValues;

    T0Layers = T0LayersPL->get(PT);
    if (T0Layers == -1 || T0Layers > y) {
        T0Layers = y;
    }
    T1Layers = T1LayersPL->get(PT);
    if (T1Layers == -1 || T1Layers > y) {
        T1Layers = y;
    }


    inputValues.resize(nrInputValues);
    outputValues.resize(T1_size); // the extra values here hold the recurrent, but can not be accessed outside the brain

    //std::cout << "AA" << std::endl;

    mutationRate_addWire = mutationRate_addWirePL->get(PT);
    mutationRate_eraseWire = mutationRate_eraseWirePL->get(PT);
    mutationRate_alterInput = mutationRate_alterInputPL->get(PT);
    mutationRate_alterOutput = mutationRate_alterOutputPL->get(PT);

    optimizeBrain = optimizeBrainPL->get(PT);

    saveBrainBehavior = saveBrainBehaviorPL->get(PT);

    outputMode = outputModePL->get(PT);

    chargeMin = chargeMinPL->get(PT);
    chargeMax = chargeMaxPL->get(PT);


    popFileColumns.clear();
    popFileColumns.push_back("wireCount");
    popFileColumns.push_back("wireCountActive");
}

// this is called when we create a new brain from scratch
WireCubeBrain::WireCubeBrain(int _nrInNodes, int _nrOutNodes, std::unordered_map<std::string, std::shared_ptr<AbstractGenome>>& _genomes,
    std::shared_ptr<ParametersTable> PT_) : WireCubeBrain(_nrInNodes, _nrOutNodes, PT_) {
}

void WireCubeBrain::update(){
    //std::cout << "in update..." << std::endl;
    if (clearBetweenUpdates) {
        std::fill(matrix_state.begin(), matrix_state.end(), 0);
    }
    //std::fill(matrix_state_next.begin(), matrix_state_next.end(), 0);

    //std::cout << "in update... 1" << std::endl;

    for (int i = 0; i < T0_size; i++) { // set the values of the cells connected to inputs and recurrent inputs
        // if i < nrInputValues, pull value from inputValues[i] (where -1 = chanrged, and 0 = not charged.
        // if i >= nrInputValues, pull from T+1 buffer (i.e., hidden in outputValues which start at nrOutputValues)
        matrix_state[T0_cells[i]] = (i < nrInputValues) ? (inputValues[i] ? -1 : 0) : (outputValues[i - nrInputValues + nrOutputValues] ? -1 : 0);
    }

    if (saveBrainBehavior == 1) { // this will only happen the first time this brain is update!
        std::string tempStr; // used if we need to save state
        tempStr = "";

        tempStr += "name,values\n";
        tempStr += "scale,\""+std::string(std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z)) + "\"\nwire,\"";
        for (auto l : matrix_wire) {
            tempStr += std::to_string(l);
        }

        std::vector<int> temp_matrix_wire_active(cellCount, 0);
        for (auto c : wireCellsActive) {
            temp_matrix_wire_active[c] = 1;
        }
        tempStr += "\"\nwire_active,\"";
        for (auto l : temp_matrix_wire_active) {
            tempStr += std::to_string(l);
        }

        tempStr += "\"\nT0_cells,\"";
        for (auto c : T0_cells) {
            tempStr += std::to_string(c)+",";
        }
        tempStr.pop_back();
        tempStr += "\"\nT1_cells,\"";
        for (auto c : T1_cells) {
            tempStr += std::to_string(c) + ",";
        }
        tempStr.pop_back();
        tempStr += "\"\nnrRecurrent," + std::to_string(nrRecurrentValues) + "\n";
        tempStr += "clearBetweenUpdates," + std::to_string(clearBetweenUpdates) + "\n";
        tempStr += "decay," + std::to_string(decayDuration) + "\n";
        tempStr += "brainUpdates," + std::to_string(stepsPerUpdate) + "\n";
        FileManager::openAndWriteToFile("config.csv", tempStr);
        saveBrainBehavior = 2; // this indicates we are on world update 0 (saveBrainBehavior - 2)
    }

        //save brain behavior
    if (saveBrainBehavior > 1) {
        std::string tempStr; // used if we need to save state
        // save update,step,T0_cells
        tempStr = "";
        for (auto c : T0_cells) {
            tempStr += std::to_string(matrix_state[c]) + ",";
        }
        tempStr.pop_back();
        FileManager::openAndWriteToFile("activations.csv", tempStr);

    }

    std::vector<double> outValuesSums(T1_size, 0.0);

    for (size_t u = 0; u < stepsPerUpdate; u++) { // for number of brain updates per world update
        ////////////////////////////////////////////////////////////////////////////////////
        // brainStateLists and fullStateBuffer were added for the WIRE_MEMORY project!
        ////////////////////////////////////////////////////////////////////////////////////
        // this is where we will save all fullStateBuffer related stuff
        // this means that the states will represent the state just before an update
        // ... so input cells will be representative of their true state.
        // This also means that we will not record the last state (which I think is Okay)
        ////////////////////////////////////////////////////////////////////////////////////
        if (saveBrainBehavior > 1) {
            std::string brainStateStr;
            for (auto& currentIndex : wireCellsActive) { // for every cell in wireCellsActive list (optimized list), determin next time step state 
                fullStateBuffer[currentIndex].NeighborStates.push_back({});
                for (auto& neighborIndex : links[currentIndex]) {
                    fullStateBuffer[currentIndex].NeighborStates.back().push_back(matrix_state[neighborIndex]);
                }
                fullStateBuffer[currentIndex].focalStates.push_back(matrix_state[currentIndex]);
                brainStateStr += std::to_string(matrix_state[currentIndex]);
            }
            brainStateLists.push_back(brainStateStr);
        }

        inputStates.push_back(inputValues);
        std::string tmpStr = "";
        for (auto v : inputValues) {
            tmpStr+=std::to_string(Bit(v));
            //std::cout << v << "    " << tmpStr << std::endl;
        }
        inputStatesStrings.push_back(tmpStr);

        //std::string brainStateStr;
        for (auto& currentIndex : wireCellsActive) { // for every cell in wireCellsActive list (optimized list), determin next time step state 
            //std::cout << "in update... wire cell with index: " << i << std::endl;

            if (matrix_state[currentIndex] == 0) { // if this cell is in the wire state...
                //std::cout << "  wire ";

                //if (saveBrainBehavior) {
                //    fullStateBuffer[currentIndex].NeighborStates.push_back({});
                //}
                size_t neighborSum = 0; // we need count the number of charged neighbors
                for (auto& neighborIndex : links[currentIndex]) {
                    neighborSum += (matrix_state[neighborIndex] == -1); // +1 for each charged neighbor
                    //std::cout << j << " " << matrix_state[j] << " " << neighborSum << " | ";
                    //if (saveBrainBehavior) {
                    //    fullStateBuffer[currentIndex].NeighborStates.back().push_back(matrix_state[neighborIndex] == -1);
                    //}
                }
                if (neighborSum > (chargeMin-1) && neighborSum < (chargeMax+1)) { // i.e., 1 or 2 charged neighbors
                    //std::cout << " -> charged!" << std::endl;
                    matrix_state_next[currentIndex] = -1; // this cell will be charge next update
                }
                else { // there is either too much or too little charge
                    //std::cout << " -> still wire!\n";
                    matrix_state_next[currentIndex] = 0; // keep current location as wire
                }
            }
            else if (matrix_state[currentIndex] == -1) { // if this cell is change
                matrix_state_next[currentIndex] = decayDuration; // this cell will be on cooldown next update
                //std::cout << " charge -> " << matrix_state_next[i] << std::endl;
            }
            else if (matrix_state[currentIndex] > 0) { // this cell is in cool down
                matrix_state_next[currentIndex] = matrix_state[currentIndex] - 1; // cooling down - if this gets to 0, then the cell is wire again, and ready to take charge next update
                //std::cout << " decay -> " << matrix_state[i] << std::endl;
            }
            else {
                std::cout << "during update, found unexpected state in matrix_state value " << matrix_state[currentIndex] << " at index " << currentIndex << ". exiting..." << std::endl;
                exit(1);
            }
            //if (saveBrainBehavior) {
            //    fullStateBuffer[currentIndex].focalStates.push_back(matrix_state_next[currentIndex]);
            //    brainStateStr += std::to_string(matrix_state_next[currentIndex]);
            //}

        }
        //if (saveBrainBehavior) {
        //    brainStateLists.push_back(brainStateStr);
        //}
        if (recordActivity) {
            lifeTimes.back()++;
        }

        matrix_state = matrix_state_next;

        if (outputMode == 1) {
            for (int index = 0; index < T1_size; index++) {
                if (matrix_wire[T1_cells[index]] == 1) { // if this location is wire (it's possible that mutation changed location, but there is still a wire here)
                    outValuesSums[index] += (matrix_state[T1_cells[index]] == -1) ? 1.0 : 0.0;
                }
            }
        }

        //lets try reseting the inputs every update. This way, they will only be on when set, and they can not be used in the computation otherwise
        for (int i = 0; i < T0_size; i++) { // set the values of the cells connected to inputs and recurrent inputs
            matrix_state[T0_cells[i]] = 0;
        }

    }
    //std::cout << "in update... done update" << std::endl;
    
    if (saveBrainBehavior > 1) { // update world update counter
        saveBrainBehavior++;
    }

    for (int index = 0; index < T1_size; index++) {
        if (outputMode == 0) { // last only
            outputValues[index] = matrix_state[T1_cells[index]] == -1;
        }
        else if (outputMode == 1) { // average over update
            outputValues[index] = (outValuesSums[index] >= ((double)stepsPerUpdate / 2.0)) ? 1 : 0;
        }
    }
}

// make a copy of the brain that called this
std::shared_ptr<AbstractBrain> WireCubeBrain::makeCopy(std::shared_ptr<ParametersTable> PT) {
    //std::cout << "in makeCopy" << std::endl;
    // You need to define this function. It needs to return a copy of the brain that called it
    auto newBrain = std::make_shared<WireCubeBrain>(nrInputValues, nrOutputValues, PT);
    newBrain->cellCount = cellCount;
    //std::cout << "   cellCount: " << cellCount << std::endl;
    //std::cout << " T0_cells.size(): " << T0_cells.size() << "   T0_size: " << T0_size << std::endl;
    //std::cout << " T1_cells.size(): " << T1_cells.size() << "   T1_size: " << T1_size << std::endl;

    newBrain->matrix_wire = matrix_wire; // size(x*y*z), matrix before updates begin
    newBrain->T0_cells = T0_cells; // indices of T0 cells (i.e., input + recurrent)
    newBrain->T1_cells = T1_cells; // indices of T1 cells (i.e., output + recurrent)

    newBrain->wireCells = wireCells; // list of indices that are wire in matrix_initial
    newBrain->links = links; // list of links (i.e., neighboors for each cell)

    newBrain->matrix_state = std::vector<int>(cellCount, 0); // we are now ready for the first update
    newBrain->matrix_state_next = std::vector<int>(cellCount, 0); // we are now ready for the first update
    //std::cout << "  done copy" << std::endl;
    return(newBrain);
}

// Make a brain like the brain that called this function, using genomes and initalizing other elements.
std::shared_ptr<AbstractBrain> WireCubeBrain::makeBrain(std::unordered_map<std::string, std::shared_ptr<AbstractGenome>> & _genomes) {
        
    // You need to define this function. It needs to return a brain, it can use the brain that called this, and _genomes
    auto newBrain = std::make_shared<WireCubeBrain>(nrInputValues, nrOutputValues, _genomes, PT);
    // create a brain and fill it in with some wire
    newBrain->matrix_wire = std::vector<int>(newBrain->cellCount, 0); // make an empty matrix
    newBrain->T0_cells.resize(newBrain->T0_size); // reserve space for inputs
    newBrain->T1_cells.resize(newBrain->T1_size); // reserve space for outputs

    size_t index;
    int intiCount = double(newBrain->cellCount) * double(newBrain->initialFillRatio);
    for (size_t i = 0; i < intiCount; i++) { // initialFillRatio * cellCount cells to 1 (wire)
        do {
            index = Random::getIndex(newBrain->cellCount); // get a random location
        } while (newBrain->matrix_wire[index] == 1); // if already wire, pick another location
        //std::cout << " ... " << index << std::endl;
        newBrain->matrix_wire[index] = 1; // set selected location to wire
        newBrain->wireCells.push_back(index); // add selected location to wireCells (list of cells that are wire)
    }

    int layerSize = newBrain->x * newBrain->z;
    // pick input and output matrix indices from wire cells
    for (size_t i = 0; i < newBrain->T0_size; i++) {
        newBrain->T0_cells[i] = Random::getIndex(layerSize * newBrain->T0Layers);
        if (newBrain->matrix_wire[newBrain->T0_cells[i]] == 0) {
            newBrain->matrix_wire[newBrain->T0_cells[i]] = 1;
            newBrain->wireCells.push_back(newBrain->T0_cells[i]);
        }
        //std::cout << "setting up T0 cell:" << i << " : " << T0_cells[i] << std::endl;
    }

    for (size_t i = 0; i < newBrain->T1_size; i++) {
        //newBrain->T1_cells[i] = newBrain->wireCells[Random::getIndex(newBrain->wireCells.size())];
        //std::cout << " BB " << tempIndex << " -- " << newBrain->T1_size << " : " << newBrain->T1_cells.size() << std::endl;
        newBrain->T1_cells[i] = Random::getInt(cellCount - (layerSize * newBrain->T1Layers), cellCount - 1);
        if (newBrain->matrix_wire[newBrain->T1_cells[i]] == 0) {
            newBrain->matrix_wire[newBrain->T1_cells[i]] = 1;
            newBrain->wireCells.push_back(newBrain->T1_cells[i]);
        }
        //std::cout << "setting up T1 cell:" << i << " : " << T1_cells[i] << std::endl;
    }

    newBrain->makeLinks();

    newBrain->matrix_state = std::vector<int>(newBrain->cellCount, 0); // we are now ready for the first update
    newBrain->matrix_state_next = std::vector<int>(newBrain->cellCount, 0); // we are now ready for the first update
    //std::cout << "done setup" << std::endl;

    return(newBrain);
}

std::string WireCubeBrain::description() {
    // returns a desription of this brain in it's current state
    return "WireCubeBrain - no description";
}

DataMap WireCubeBrain::getStats(std::string& prefix) {
    // return a vector of DataMap of stats from this brain, this is called just after the brain is constructed
    // values in this datamap are added to the datamap on the organism that "owns" this brain
    // all data names must have prefix prepended (i.e. connections would be prefix + "connections"
    
    DataMap dataMap;
    // datamap example:
    //dataMap.append(prefix + "someStatName", someStat);
    dataMap.append(prefix + "wireCount", (int)wireCells.size());
    dataMap.append(prefix + "wireCountActive", (int)wireCellsActive.size());
    return dataMap;
}

std::string WireCubeBrain::getType() {
    // return the type of this brain
    return "WireCubeBrain";
}

void WireCubeBrain::setInput(const int& inputAddress, const double& value) {
    if (inputAddress < nrInputValues) {
        inputValues[inputAddress] = value;
    }
    else {
        std::cout << "in Brain::setInput() : Writing to invalid input ("
            << inputAddress << ") - this brain needs more inputs!\nExiting"
            << std::endl;
        exit(1);
    }
}

double WireCubeBrain::readInput(const int& inputAddress) {
    if (inputAddress < nrInputValues) {
        return inputValues[inputAddress];
    }
    else {
        std::cout << "in Brain::readInput() : Reading from invalid input ("
            << inputAddress << ") - this brain needs more inputs!\nExiting"
            << std::endl;
        exit(1);
    }
}

void WireCubeBrain::setOutput(const int& outputAddress, const double& value) {
    if (outputAddress < nrOutputValues) {
        outputValues[outputAddress] = value;
    }
    else {
        std::cout << "in Brain::setOutput() : Writing to invalid output ("
            << outputAddress << ") - this brain needs more outputs!\nExiting"
            << std::endl;
        exit(1);
    }
}

double WireCubeBrain::readOutput(const int& outputAddress) {
    if (outputAddress < nrOutputValues) {
        return outputValues[outputAddress];
    }
    else {
        std::cout << "in Brain::readOutput() : Reading from invalid output ("
            << outputAddress << ") - this brain needs more outputs!\nExiting"
            << std::endl;
        exit(1);
    }
}

void WireCubeBrain::resetOutputs() {
    for (int i = 0; i < T1_size; i++) { // clear all values (includes recurrent)
        outputValues[i] = 0.0;
    }
}

void WireCubeBrain::resetInputs() {
    for (int i = 0; i < nrInputValues; i++) {
        inputValues[i] = 0.0;
    }
}

void WireCubeBrain::resetBrain() {
    resetInputs();
    resetOutputs();

    // we need to clearn this incase clear between updates is off!
    std::fill(matrix_state.begin(), matrix_state.end(), 0);

    if (recordActivity) {
        if (lifeTimes.back() != 0) {
            lifeTimes.push_back(0);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// these functions need to be filled in if genomes are being used in this brain
///////////////////////////////////////////////////////////////////////////////////////////

// return a set of namespaces, MABE will insure that genomes with these names are created
// on organisms with these brains.
std::unordered_set<std::string> WireCubeBrain::requiredGenomes() {
    // note namespace must end with ::, also, if you wish to use default values, then 
    //    return {"root::"};

    // this brain has no king, this brain needs no king! er, genome...
    return {};
}

void WireCubeBrain::initializeGenomes(std::unordered_map<std::string, std::shared_ptr<AbstractGenome>> & _genomes) {
    // do nothing by default... since this is a direct encoded brain, then no action is needed.
    // This can be used to randomize the genome and/or insert start codons
    // genomes will be found in _genomes[name] where name is the string used in required genomes
    // example: _genomes["root::"]

}

///////////////////////////////////////////////////////////////////////////////////////////
// these functions need to be filled in if this brain is direct encoded (in part or whole)
///////////////////////////////////////////////////////////////////////////////////////////

// Make a brain like the brain that called this function, using genomes and
// inheriting other elements from parent.
std::shared_ptr<AbstractBrain> WireCubeBrain::makeBrainFrom(
    std::shared_ptr<AbstractBrain> parent,
    std::unordered_map<std::string, 
    std::shared_ptr<AbstractGenome>> & _genomes) {
   
    // in the default case, we assume geneticly encoded brains, so this just calls
    // the no parent version (i.e. makeBrain which builds from genomes)

    // since this is direct encoded, we call make copy - MABE will call mutate() later
    return makeCopy(PT);
}

// see makeBrainFrom, same thing, but for more then one parent
std::shared_ptr<AbstractBrain> WireCubeBrain::makeBrainFromMany(
    std::vector<std::shared_ptr<AbstractBrain>> parents,
    std::unordered_map<std::string,
    std::shared_ptr<AbstractGenome>> & _genomes) {

    return makeCopy(PT);
}

// if index is in T0_cells or T1_cells return false
bool WireCubeBrain::indexNotLinked(int index) {
    bool clear = true;
    for (auto c : T0_cells) {
        if (index == c) {
            clear = false;
        }
    }
    for (auto c : T1_cells) {
        if (index == c) {
            clear = false;
        }
    }
    return clear;
}


// if index is in T0_cells or T1_cells remove from matrix_wire (set to 0) and remove from wireCells
void WireCubeBrain::clearIfNotLinked(int index) {
    if (indexNotLinked(index)) {
        //std::cout << "cleared: " << index << std::endl;
        auto iter = std::find(wireCells.begin(), wireCells.end(), index);
        wireCells[iter - wireCells.begin()] = wireCells.back();
        wireCells.pop_back();
        matrix_wire[index] = 0;
    }
    //else {
    //    std::cout << "XXXXXXX: " << index << std::endl;
    //}
}

// apply direct mutations to this brain
void WireCubeBrain::mutate() {
    //std::cout << "in mutate..." << std::endl;
    //report();
    
    // mutate matrix
    

    bool showmutations = false;
    if (showmutations) {
        std::cout << "--------------------------------------------------------------------" << std::endl;
        report();
    }

    // erase some
    size_t numMutations = Random::getBinomial(wireCells.size(), mutationRate_eraseWire);
    size_t index;
    size_t cellIndex;
    for (size_t i = 0; i < numMutations; i++) {
        cellIndex = Random::getIndex(wireCells.size()); // select a cell that is wire
        if (matrix_wire[wireCells[cellIndex]] == 0) {
            std::cout << "matrix_wire[wireCells[" << cellIndex << "]] should be 1" << std::endl;
            std::cout << "matrix_wire[" << wireCells[cellIndex] << "] should be 1" << std::endl;
            std::cout << "found problem!" << std::endl;
            exit(1);
        }
        if (showmutations) {
            std::cout << "mutation attempting to remove " << wireCells[cellIndex] << std::endl;
        }
        clearIfNotLinked(wireCells[cellIndex]);
    }

    if (showmutations) {
        std::cout << "removed wire("<< numMutations <<"):" << std::endl;
        report();
    }

    //add some
    numMutations = Random::getBinomial(cellCount - wireCells.size(), mutationRate_addWire);
    for (size_t i = 0; i < numMutations; i++) { // flip some bits in the matrix
        do{ // select a cell that is empty (not wire)
            index = Random::getIndex(cellCount);
        } while (std::find(wireCells.begin(), wireCells.end(), index) != wireCells.end());
        if (matrix_wire[index] == 1) {
            std::cout << "matrix_wire[" << index << "] should be 0" << std::endl;
            std::cout << "found problem!" << std::endl;
            exit(1);
        }
        if (showmutations) {
            std::cout << "mutation has added " << index << std::endl;
        }
        matrix_wire[index] = 1;
        wireCells.push_back(index);
    }

    if (showmutations) {
        std::cout << "added wire(" << numMutations << "):" << std::endl;
        report();
    }
    //std::cout << " done matrix mutations" << std::endl;

    int layerSize = x * z;
    // mutate T0
    numMutations = Random::getBinomial(T0_size, mutationRate_alterInput);
    for (size_t i = 0; i < numMutations; i++) { // rewire some T0 indices
        //std::cout << " i: " << i << " T0_cells.size(): " << T0_cells.size() << "   T0_size: " << T0_size << std::endl;
        //T0_cells[Random::getIndex(T0_size)] = wireCells[Random::getIndex(wireCells.size())];
        cellIndex = Random::getIndex(T0_size); // pick an input cell to change
        size_t oldIndex = T0_cells[cellIndex]; // save the wire cube inde input is currently wired to
        do { // pick a new differnt input cell
            if (showmutations) {
                std::cout << cellIndex << "(" << T0_cells[cellIndex] << "/" << oldIndex << ") : " << layerSize << " * " << T0Layers << " = " << layerSize * T0Layers << std::endl;
            }
            T0_cells[cellIndex] = Random::getIndex(layerSize * T0Layers);
        } while (T0_cells[cellIndex] == oldIndex);
        if (matrix_wire[T0_cells[cellIndex]] == 0) {
            matrix_wire[T0_cells[cellIndex]] = 1;
            wireCells.push_back(T0_cells[cellIndex]);
        }
        clearIfNotLinked(oldIndex);
    }
    //std::cout << " done T0_cells mutate" << std::endl;
    if (showmutations) {
        std::cout << "T0 updated(" << numMutations << "):" << std::endl;
        report();
    }

    // mutate T1
    numMutations = Random::getBinomial(T1_size, mutationRate_alterOutput);
    for (size_t i = 0; i < numMutations; i++) { // rewire some T1 indices
        //std::cout << " i: " << i << " T1_cells.size(): " <<  T1_cells.size() << "   T1_size: " << T1_size << std::endl;
        //T1_cells[Random::getIndex(T1_size)] = wireCells[Random::getIndex(wireCells.size())];
        cellIndex = Random::getIndex(T1_size); // pick an output to change
        size_t oldIndex = T1_cells[cellIndex]; // save the wire cube index output is currently wired to
        do { // pick a different new output cell
            if (showmutations) {
                std::cout << cellIndex << "(" << T1_cells[cellIndex] << "/" << oldIndex << ") : (" << cellCount - (layerSize * T1Layers) << "," << cellCount - 1 << ")" << std::endl;
            }
            T1_cells[cellIndex] = Random::getInt(cellCount - (layerSize * T1Layers), cellCount - 1);
        } while (T1_cells[cellIndex] == oldIndex);
        if (matrix_wire[T1_cells[cellIndex]] == 0) { // if it's not wire, make it wire
            matrix_wire[T1_cells[cellIndex]] = 1;
            wireCells.push_back(T1_cells[cellIndex]);
        }
        clearIfNotLinked(oldIndex);
    }
    if (showmutations) {
        std::cout << "T1 updated(" << numMutations << "):" << std::endl;
        report();
    }

    //std::cout << " done T1_cells mutate" << std::endl;

    matrix_state = std::vector<int>(cellCount, 0); // we are now ready for the first update
    matrix_state_next = std::vector<int>(cellCount, 0); // we are now ready for the first update

    //std::cout << " done mutate" << std::endl;

    makeLinks();

    //std::cout << " done update links" << std::endl;

    //report
    if (showmutations) {
        std::cout << "DONE:" << std::endl;
        report();
    }
    // do nothing by default...
    // if this is a direct encoded brain, then this function needs to be filled in to allow for mutations.
}

// convert a brain into data map with data that can be saved to file so this brain can be reconstructed
// 'name' here contains the prefix that must be so that deserialize can identify relavent data
DataMap WireCubeBrain::serialize(std::string & name) {
    DataMap dataMap;
    std::string outStr = "";

    for (auto l : matrix_wire) {
        outStr += std::to_string(l);
    }
    outStr += "|";
    for (auto c : T0_cells) {
        outStr += std::to_string(c) + "^";
    }
    outStr.pop_back();
    outStr += "|";
    for (auto c : T1_cells) {
        outStr += std::to_string(c) + "^";
    }
    outStr.pop_back();

    //std::cout << outStr << std::endl;
    dataMap.set(name + "brainData", outStr);
    return dataMap;
}

// given an unordered_map<string, string> of org data and PT, load data into this brain
// 'name' here contains the prefix that was used when data was being saved
void WireCubeBrain::deserialize(std::shared_ptr<ParametersTable> PT,
    std::unordered_map<std::string, std::string> & orgData,
    std::string & name) {
    // note that the process by which deserialization (including text formatting) depends on
    // the corisponding serialize process    


    std::cout << "In WireCubeBrain::deserialize" << std::endl;

    std::string data = orgData[name + "brainData"];
    
    std::vector<std::string> BrainStr; // a list with all the brains parts
    convertCSVListToVector(data, BrainStr, '|'); // first division of data, matrix,T0_cells, and T1_cells

    matrix_wire = std::vector<int>(cellCount, 0); // make an empty matrix
    wireCells.clear();
    T0_cells.resize(T0_size); // reserve space for inputs
    T1_cells.resize(T1_size); // reserve space for outputs

    for (size_t i = 0; i < cellCount; i++) {
        if (BrainStr[0][i] == '1') {
            matrix_wire[i] = BrainStr[0][i] == '1';
            wireCells.push_back(i);
        }
    }
    convertCSVListToVector(BrainStr[1], T0_cells, '^');
    convertCSVListToVector(BrainStr[2], T1_cells, '^');

    makeLinks();

    matrix_state = std::vector<int>(cellCount, 0); // we are now ready for the first update
    matrix_state_next = std::vector<int>(cellCount, 0); // we are now ready for the first update

    report();

    auto test_serialize = serialize(name);
    test_serialize.openAndWriteToFile("WireSerializeTest.csv");
    std::cout << "saved a copy to WireSerializeTest.csv" << std::endl;
}
