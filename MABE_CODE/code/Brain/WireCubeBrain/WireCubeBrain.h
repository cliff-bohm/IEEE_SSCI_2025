//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#pragma once					// directive to insure that this .h file is only included one time

// AbstractBrain defines all the basic function templates for brains
#include <Brain/AbstractBrain.h>

// If your brain is (or maybe) constructed using a genome, you must include AbstractGenome.h
#include <Genome/AbstractGenome.h>

class WireCubeBrain : public AbstractBrain {

    static std::shared_ptr<ParameterLink<int>> xPL;
    int x;
    static std::shared_ptr<ParameterLink<int>> yPL;
    int y;
    static std::shared_ptr<ParameterLink<int>> zPL;
    int z;
    static std::shared_ptr<ParameterLink<int>> decayDurationPL;
    int decayDuration;
    static std::shared_ptr<ParameterLink<int>> stepsPerUpdatePL;
    int stepsPerUpdate;
    static std::shared_ptr<ParameterLink<bool>> clearBetweenUpdatesPL;
    bool clearBetweenUpdates;
    static std::shared_ptr<ParameterLink<double>> initialFillRatioPL;
    double initialFillRatio;

    static std::shared_ptr<ParameterLink<double>> mutationRate_addWirePL;
    double mutationRate_addWire;
    static std::shared_ptr<ParameterLink<double>> mutationRate_eraseWirePL;
    double mutationRate_eraseWire;
    static std::shared_ptr<ParameterLink<double>> mutationRate_alterInputPL;
    double mutationRate_alterInput;
    static std::shared_ptr<ParameterLink<double>> mutationRate_alterOutputPL;
    double mutationRate_alterOutput;

    

    static std::shared_ptr<ParameterLink<bool>> optimizeBrainPL;
    bool optimizeBrain;

    static std::shared_ptr<ParameterLink<int>> saveBrainBehaviorPL;
    int saveBrainBehavior;

    static std::shared_ptr<ParameterLink<int>> outputModePL;
    int outputMode;

    static std::shared_ptr<ParameterLink<int>> chargeMinPL;
    int chargeMin;
    static std::shared_ptr<ParameterLink<int>> chargeMaxPL;
    int chargeMax;

    static std::shared_ptr<ParameterLink<int>> T0LayersPL;
    int T0Layers;
    static std::shared_ptr<ParameterLink<int>> T1LayersPL;
    int T1Layers;

    size_t T0_size, T1_size; // size of T and T+1 buffers

    //static std::shared_ptr<ParameterLink<bool>> wormholeRatePL;
    //static std::shared_ptr<ParameterLink<int>> wiregenesWormholesBidirectionalPL;

    static std::shared_ptr<ParameterLink<int>> nrRecurrentValuesPL;
    int nrRecurrentValues;

    size_t cellCount;

    std::vector<int> matrix_wire; // size(x*y*z), structure of brain, only used durring setup
    std::vector<int> T0_cells; // indices of T0 cells (i.e., input + recurrent)
    std::vector<int> T1_cells; // indices of T1 cells (i.e., output + recurrent)

    std::vector<int> wireCells; // list of indices that are wire in matrix_initial
    std::vector<std::vector<int>> links; // list of links (i.e., neighboors for each cell)

    std::vector<int> matrix_state; //                        <--<--<--<--<
    std::vector<int> matrix_state_next; // matrix we update into -->-->--^

public:
    std::vector<int> wireCellsActive; // list of indices that are wire in matrix_initial

    WireCubeBrain() = delete;


    class WireStateBuffer { // this will store the state of one focal node and the states of the neigbors from the prior time step
    public:
        std::vector<int> focalStates; // post state for every wire
        std::vector<std::vector<int>> NeighborStates; // pre states for every neighbor of this focal, each list is one brain update
        int ID = -1; //id for this node (i.e., locator in brain)
        std::vector<int> neighborIDs;
    };

    std::vector<WireStateBuffer> fullStateBuffer; // states buffers for all wire in brain
    std::vector<std::vector<double>> inputStates; // inputValues (will remain constant over all brain updates of a given world update)
    std::vector<std::string> inputStatesStrings; // inputValues (will remain constant over all brain updates of a given world update)
    std::vector<std::string> brainStateLists;


    WireCubeBrain(int ins, int outs, std::shared_ptr<ParametersTable> PT);
    WireCubeBrain(int _nrInNodes, int _nrOutNodes,
        std::unordered_map<std::string, std::shared_ptr<AbstractGenome>>& _genomes,
        std::shared_ptr<ParametersTable> PT_ = nullptr);

    virtual ~WireCubeBrain() = default;

    void makeLinks();
    void reportFull();
    void report();
    void showBrainActivityState();
    bool indexNotLinked(int index);
    void clearIfNotLinked(int index);

    virtual void update();

    // make a copy of the brain that called this
    virtual std::shared_ptr<AbstractBrain> makeCopy(std::shared_ptr<ParametersTable> PT);

    // Make a brain like the brain that called this function, using genomes and initalizing other elements.
    virtual std::shared_ptr<AbstractBrain> makeBrain(std::unordered_map<std::string, std::shared_ptr<AbstractGenome>>& _genomes);

    virtual std::string description(); // returns a desription of this brain in it's current state

    virtual DataMap getStats(std::string& prefix); // return a vector of DataMap of stats from this brain

    virtual std::string getType(); // return the type of this brain

    virtual void setInput(const int& inputAddress, const double& value);

    virtual double readInput(const int& inputAddress);

    virtual void setOutput(const int& outputAddress, const double& value);

    virtual double readOutput(const int& outputAddress);

    virtual void resetOutputs();

    virtual void resetInputs();

    virtual void resetBrain();

    ///////////////////////////////////////////////////////////////////////////////////////////
    // these functions need to be filled in if genomes are being used in this brain
    ///////////////////////////////////////////////////////////////////////////////////////////

    virtual std::unordered_set<std::string> requiredGenomes(); // does this brain use any genomes
    
    // initializeGenomes can be used to randomize the genome and/or insert start codons
    virtual void initializeGenomes(std::unordered_map<std::string, std::shared_ptr<AbstractGenome>>& _genomes);

    ///////////////////////////////////////////////////////////////////////////////////////////
    // these functions need to be filled in if this brain is direct encoded (in part or whole)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // Make a brain like the brain that called this function, using genomes and
    // inheriting other elements from parent.
    // in the default case, we assume geneticly encoded brains, so this just calls
    // the no parent version (i.e. makeBrain which builds from genomes)
    virtual std::shared_ptr<AbstractBrain> makeBrainFrom(
        std::shared_ptr<AbstractBrain> parent,
        std::unordered_map<std::string,
        std::shared_ptr<AbstractGenome>>& _genomes);

    // see makeBrainFrom, this can take more then one parent
    virtual std::shared_ptr<AbstractBrain> makeBrainFromMany(
        std::vector<std::shared_ptr<AbstractBrain>> parents,
        std::unordered_map<std::string,
        std::shared_ptr<AbstractGenome>>&_genomes);

    // apply direct mutations to this brain
    virtual void mutate();

    // convert a brain into data map with data that can be saved to file so this brain can be reconstructed
    // 'name' here contains the prefix that must be so that deserialize can identify relavent data
    virtual DataMap serialize(std::string& name);

    // given an unordered_map<string, string> of org data and PT, load data into this brain
    // 'name' here contains the prefix that was used when data was being saved
    virtual void deserialize(std::shared_ptr<ParametersTable> PT,
        std::unordered_map<std::string, std::string>& orgData,
        std::string& name);

};

inline std::shared_ptr<AbstractBrain> WireCubeBrain_brainFactory(int ins, int outs, std::shared_ptr<ParametersTable> PT) {
    return std::make_shared<WireCubeBrain>(ins, outs, PT);
}
