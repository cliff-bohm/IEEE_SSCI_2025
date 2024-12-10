//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

// Evaluates agents on how many '1's they can output. This is a purely fixed
// task
// that requires to reactivity to stimuli.
// Each correct '1' confers 1.0 point to score, or the decimal output determined
// by 'mode'.

#include "NBackWorld.h"
#include "../../Brain/WireCubeBrain/WireCubeBrain.h"


std::shared_ptr<ParameterLink<int>> NBackWorld::testMutantsPL = Parameters::register_parameter("WORLD_NBACK-testMutants", 0, "if > 0, this number of mutants of each agent will be tested");

std::shared_ptr<ParameterLink<int>> NBackWorld::evaluationsPerGenerationPL =
Parameters::register_parameter("WORLD_NBACK-evaluationsPerGeneration", 10, "Number of times to evaluate each agent per generation");

std::shared_ptr<ParameterLink<int>> NBackWorld::testsPerEvaluationPL =
Parameters::register_parameter("WORLD_NBACK-testsPerEvaluation", 10, "Number of times to test each agent per evaluation");

std::shared_ptr<ParameterLink<std::string>> NBackWorld::NsListsPL =
Parameters::register_parameter("WORLD_NBACK-NsList", (std::string)"1,2,3:100|2,3,4:-1", "comma seperated list of n values followed by ':' and a time\n"
	"more then one list can be defined seperated by '|'. The last list time must be -1 (i.e. forever)\n"
	"eg: 1,2,3:100|2,3,4:-1");

std::shared_ptr<ParameterLink<int>> NBackWorld::scoreMultPL =
Parameters::register_parameter("WORLD_NBACK-scoreMult", 1, "score multiplier");

std::shared_ptr<ParameterLink<int>> NBackWorld::RMultPL =
Parameters::register_parameter("WORLD_NBACK-RMult", 1, "score R multiplier");

std::shared_ptr<ParameterLink<int>> NBackWorld::delayOutputEvalPL =
Parameters::register_parameter("WORLD_NBACK-delayOutputEval", 0, "generation delay for ouput evalutation");

std::shared_ptr<ParameterLink<bool>> NBackWorld::tritInputsPL =
Parameters::register_parameter("WORLD_NBACK-tritInputs", false, "if false (defaut) then inputs to brain are 0 or 1. If true, inputs are -1,0,1");

std::shared_ptr<ParameterLink<std::string>> NBackWorld::groupNamePL =
Parameters::register_parameter("WORLD_NBACK-groupNameSpace", (std::string) "root::", "namespace of group to be evaluated");
std::shared_ptr<ParameterLink<std::string>> NBackWorld::brainNamePL =
Parameters::register_parameter("WORLD_NBACK-brainNameSpace", (std::string) "root::", "namespace for parameters used to define brain");

std::shared_ptr<ParameterLink<bool>> NBackWorld::saveFragOverTimePL =
Parameters::register_parameter("WORLD_NBACK_ANALYZE-saveFragOverTime", false,
	"");
std::shared_ptr<ParameterLink<bool>> NBackWorld::saveBrainStructureAndConnectomePL =
Parameters::register_parameter("WORLD_NBACK_ANALYZE-saveBrainStructureAndConnectome", true,
	"");
std::shared_ptr<ParameterLink<bool>> NBackWorld::saveStateToStatePL =
Parameters::register_parameter("WORLD_NBACK_ANALYZE-saveStateToState", true,
	"");
std::shared_ptr<ParameterLink<bool>> NBackWorld::save_R_FragMatrixPL =
Parameters::register_parameter("WORLD_NBACK_ANALYZE-save_R_FragMatrix", false,
	"");
std::shared_ptr<ParameterLink<bool>> NBackWorld::saveFlowMatrixPL =
Parameters::register_parameter("WORLD_NBACK_ANALYZE-saveFlowMatrix", false,
	"");
std::shared_ptr<ParameterLink<bool>> NBackWorld::saveStatesPL =
Parameters::register_parameter("WORLD_NBACK_ANALYZE-saveStates", true,
	"");


#include "../../Utilities/PowerSet.h"

NBackWorld::NBackWorld(std::shared_ptr<ParametersTable> PT_) : AbstractWorld(PT_) {

	saveFragOverTime = saveFragOverTimePL->get(PT);
	saveBrainStructureAndConnectome = saveBrainStructureAndConnectomePL->get(PT);
	saveStateToState = saveStateToStatePL->get(PT);
	save_R_FragMatrix = save_R_FragMatrixPL->get(PT);
	saveFlowMatrix = saveFlowMatrixPL->get(PT);
	saveStates = saveStatesPL->get(PT);

	std::vector<std::string> NListsBreakDown1; // used to parse nLists
	std::vector<std::string> NListsBreakDown2; // used to parse nLists

	convertCSVListToVector(NsListsPL->get(PT), NListsBreakDown1, '|'); // lists (i.e. Ns to score + time) are sperated by '|'

	int temp = 0;

	for (auto elem : NListsBreakDown1) {
		convertCSVListToVector(elem, NListsBreakDown2, ':'); // list of Ns is sperated from time with a ':'
		convertString(NListsBreakDown2[1], temp); // get the time for this list
		if (NListSwitchTimes.size() == 0) { // if this is the first list, then put this time on NListSwitchTimes
			NListSwitchTimes.push_back(temp);
		}
		else if (temp > 0) { // else, if it's not -1 (i.e. last list), put this time + sum of previous times
			NListSwitchTimes.push_back(NListSwitchTimes.back() + temp);
		}
		else { // else it's -1, this is the last list (note, if the first list has time -1, that is handled by the if)
			NListSwitchTimes.push_back(-1);
		}
		NListLists.push_back({}); // add a blank list so we have a container to fill
		convertCSVListToVector(NListsBreakDown2[0], NListLists.back(), ','); // fill the container we just added to NListLists
	}
	
	int nextOut = 0;
	
	std::cout << "testing Lists will change on updates:";
	for (auto elem : NListSwitchTimes) {
		std::cout << "  " << elem;
	}
	std::cout << std::endl;

	std::cout << "testing Lists:\n";
	for (auto elem : NListLists) {
		for (auto elem2 : elem) {
			std::cout << "  " << elem2;
			largestN = std::max(largestN, elem2);
			if (!N2OutMap.count(elem2)) {
				N2OutMap[elem2] = nextOut++;
			}
		}
		std::cout << std::endl;
	}

	std::cout << "  largest N found was : " << largestN << ". Brains will be run for this number of world steps before testing begins." << std::endl;

	// now get currentLargestN

	for (auto elem : NListLists[0]) {
		currentLargestN = std::max(currentLargestN, elem);
	}

	evaluationsPerGeneration = evaluationsPerGenerationPL->get(PT); // each agent sees this number of inputs (+largest N) and is scored this number of times each evaluation
	testsPerEvaluation = testsPerEvaluationPL->get(PT); // each agent is reset and evaluated this number of times

	std::cout << "output map:\n";
	for (auto elem : N2OutMap) {
		std::cout << "  N: " << elem.first << " <- output: " << elem.second << std::endl;
	}
	std::cout << "brains will have 1 input and " << N2OutMap.size() << " outputs." << std::endl;

	groupName = groupNamePL->get(PT);
	brainName = brainNamePL->get(PT);



	std::vector<std::string> temp_analyzeWhatStr;


	testMutants = testMutantsPL->get(PT);

	tritInputs = tritInputsPL->get(PT);

	// columns to be added to ave file
	popFileColumns.clear();
	popFileColumns.push_back("score");
	//////popFileColumns.push_back("R");
	for (auto elem : N2OutMap) {
		popFileColumns.push_back("nBack_" + std::to_string(elem.first));
		std::cout << "adding: " << "nBack_" + std::to_string(elem.first) << std::endl;
	}
}

void NBackWorld::evaluate(std::map<std::string, std::shared_ptr<Group>> &groups, int analyze, int visualize, int debug) {

	// check to see if we need to advance to next NsList
	if (Global::update >= NListSwitchTimes[currentNList] && NListSwitchTimes[currentNList] != -1) {
		currentNList += 1;
		//std::cout << "advancing to next list... " << currentNList << std::endl;
		currentLargestN = 0;
		for (auto elem : NListLists[currentNList]) {
			currentLargestN = std::max(currentLargestN, elem);
		}
	}

	int popSize = groups[groupNamePL->get(PT)]->population.size();
	for (int i = 0; i < popSize; i++) {
		// eval this agent
		evaluateSolo(groups[groupNamePL->get(PT)]->population[i], analyze, visualize, debug);
		// now lets test some god damn dirty mutants!

		if (testMutants > 0) {
			std::vector<double> mutantScores;
			double mutantScoreSum = 0;
			for (int j = 0; j < testMutants; j++) {
				auto mutantOffspring = groups[groupNamePL->get(PT)]->population[i]->makeMutatedOffspringFrom(groups[groupNamePL->get(PT)]->population[i]);
				evaluateSolo(mutantOffspring, 0, 0, 0);
				auto s = mutantOffspring->dataMap.getAverage("score");
				mutantScores.push_back(s);
				mutantScoreSum += s;
			}
			std::cout << "score: " << groups[groupNamePL->get(PT)]->population[i]->dataMap.getAverage("score") << "  mutantAveScore(" << testMutants << "): " << mutantScoreSum / testMutants << std::endl;
			//std::ofstream mutantScoreFile;
			//mutantScoreFile.open("mutantScoreFile.csv",
			//	std::ios::out |
			//	std::ios::app);
			//mutantScoreFile << groups[groupNamePL->get(PT)]->population[i]->dataMap.getAverage("score") << "," << mutantScoreSum / testMutants << std::endl;
			//mutantScoreFile.close();
			FileManager::writeToFile("mutantScoreFile.txt", std::to_string(groups[groupNamePL->get(PT)]->population[i]->dataMap.getAverage("score")) + "," + std::to_string(mutantScoreSum / testMutants));


		}

	}
	
	
	//if (analyze) {
	//	groups[groupNamePL->get(PT)]->archive();
	//}
}

void NBackWorld::evaluateSolo(std::shared_ptr<Organism> org, int analyze, int visualize, int debug) {
	

	auto brain = org->brains[brainName];
	if (visualize) {
		brain->setRecordActivity(true);
	}

	double score = 0.0;
	std::vector<int> tallies(N2OutMap.size(), 0); // how many times did brain get each N in current list correct?
	std::vector<int> inputList(testsPerEvaluation + currentLargestN, 0);

	TS::intTimeSeries worldStates;
	//std::cout << "currentLargestN: " << currentLargestN << std::endl;
	for (int r = 0; r < evaluationsPerGeneration; r++) {
		brain->resetBrain();

		int t = 0;


		for (int t = 0; t < inputList.size(); t++) {

			int lowBound = 0;
			if (tritInputs) {
				lowBound = -1;
			}
			inputList[t] = Random::getInt(lowBound,1);
			brain->setInput(0, inputList[t]);

			brain->update();

			// collect score and world data but only once we have reached currentLargestN
			if (t >= currentLargestN) {
				// add space in world data vector
				worldStates.push_back({});
				//for (int i = currentLargestN; i >= 0; i--) {
				for (int i = 0; i <= currentLargestN; i++) { // this will by n-0 first rather than last
					worldStates.back().push_back(inputList[t - i]); // get all history, not just outputs
				}
				for (auto elem : NListLists[currentNList]) {
					//worldStates.back().push_back(inputList[t - elem + 1]);
					if (Global::update >= delayOutputEvalPL->get(PT) && // if update is greater than delay time
						Trit(brain->readOutput(N2OutMap[elem])) == inputList[t - elem]) { // if output is correct 
						score += 1; // add 1 to score
						tallies[N2OutMap[elem]] += 1; // add 1 to correct outputs for this N
					}
				}
			}
		}
	}
	

	org->dataMap.append("score", (score*scoreMultPL->get(PT)) / (evaluationsPerGeneration*testsPerEvaluation*NListLists[currentNList].size()));
	// score is divided by number of evals * number of tests * number of N's in current list

	for (auto elem : N2OutMap) {
		org->dataMap.append("nBack_" + std::to_string(elem.first), (double)tallies[elem.second] / (double)(evaluationsPerGeneration*testsPerEvaluation));
		// since this is for one N at a time, it's just divided by number of evals * number of tests
	}

	if (visualize) {
		std::cout << "organism with ID " << org->ID << " scored " << org->dataMap.getAverage("score") << std::endl;
	}
	
	
	if (0) {
		std::cout << "Generating State 2 State" << std::endl;

		// cast the brain to a wire cube brain so we can get at the guts!
		auto WireBrain = std::dynamic_pointer_cast<WireCubeBrain>(brain);

		std::unordered_map<string, string> symbolMap; // for each rawState, it's mapped Name
		std::unordered_map<string, int> stateTransCounts; // for each transition, how often does it occur?
		std::unordered_map<string, std::vector<int>> stateTransWhen;// for each transition, when does it occur?

		std::string tempList;
		int nextSymbolCode = 0;
		std::string lastSymbol;

		int when = 0;
		for (auto& state : WireBrain->brainStateLists) {
			if (symbolMap.find(state) == symbolMap.end()) {
				// get new symbol
				int n = nextSymbolCode++;
				std::string symbol;
				while (n > 0) {
					int remainder = n % 26;
					symbol = char('A' + remainder) + symbol;
					n /= 26;
				}
				if (symbol.empty()) symbol = "AAA";
				while (symbol.size() < 3) {
					symbol = "A" + symbol;
				}
				// get new symbol
				symbolMap[state] = symbol;
			}

			if (when > 0) {
				std::string transition = lastSymbol + "->" + symbolMap[state];
				stateTransCounts[transition]++;
				if (stateTransWhen.find(transition) == stateTransWhen.end()) {
					stateTransWhen[transition] = std::vector<int>();
				}
				stateTransWhen[transition].push_back(when);

			}
			lastSymbol = symbolMap[state];
			when += 1;

		}

		for (auto& t : stateTransCounts) {
			double hue = 0;
			for (auto& v : stateTransWhen[t.first]) {
				hue += v;
			}
			std::cout << t.first << " [color=\"" << hue / double(stateTransWhen[t.first].size()) / double(when * 2) << ",1,.8\",penwidth=" << 5.0 + double(stateTransCounts[t.first]) / 5.0 << "]" << std::endl;
		}
	}

  ////////////////////////////////////////////////////////
  // THIS IS NEW CODE TO SUPPORT WIRE BRAIN MEMORY PROJECT
  ////////////////////////////////////////////////////////

	if (visualize) {
		auto WireBrain = std::dynamic_pointer_cast<WireCubeBrain>(brain);
		std::cout << "score: " << score << std::endl;
		std::cout << "brainType: " << WireBrain->getType() << std::endl;
		std::cout << "stateCount = " << WireBrain->brainStateLists.size() << std::endl;;

		std::cout << "WireBrain->lifeTimes count: " << WireBrain->lifeTimes.size() << std::endl;
		for (auto lt : WireBrain->lifeTimes) {
			std::cout << lt << "  ";
		}
		std::cout << std::endl;
		auto lifeTimes = WireBrain->lifeTimes;


		// inflate worldStates so that we have 12 states per world update (the number of brain updates)
		// note that we only start collecting work state once we have currentLargestN inputs, so we won't need to trim this
		TS::intTimeSeries inflatedWorldStates;
		for (auto& WS : worldStates) {
			for (int i = 0; i < 12; i++) {
				inflatedWorldStates.push_back(WS);
			}
		}

		// for each cell that is wire [outer vector]
		// for each update [middle vector]
		// have a list with cell value
		std::vector<std::vector<std::vector<int>>> cellStates;

		// for each cell that is wire [outer vector]
		// for each update [middle vector]
		// have a list with prior neighbor values (up to 27)
		std::vector<std::vector<std::vector<int>>> neighborStates;

		for (auto& i : WireBrain->wireCellsActive) {
			std::cout << " running cell " << i << " full" << std::endl;

			// save id,neihborIDs (one list, first element is focal)
			std::string neihborIDstrings = "";
			for (auto& neighborID : WireBrain->fullStateBuffer[i].neighborIDs) {
				neihborIDstrings += std::to_string(neighborID) + ",";
			}
			neihborIDstrings.pop_back();
			FileManager::openAndWriteToFile("vis_output/neighbor_ids.txt", std::to_string(i) + ":" + neihborIDstrings);


			std::vector<int> decomposedLifetimes = lifeTimes;
			// lifetimes are not 1/12th as long!
			for (auto& l : decomposedLifetimes) {
				l = int(l / 12);
			}

			cellStates.push_back({});
			neighborStates.push_back({});
			for (int j = 0; j < WireBrain->fullStateBuffer[i].focalStates.size(); j++) {

				/////// TRYING OUT, for prediction, we only want to try to predict focal state on vs off (so all values not == -1 are considered off)
				cellStates.back().push_back({ WireBrain->fullStateBuffer[i].focalStates[j] == -1 });
				neighborStates.back().push_back({});
				for (auto& neighborID : WireBrain->fullStateBuffer[i].neighborIDs) {

				/////// TRYING OUT, for prediction, we only want to try to predict focal state on vs off (so all values not == -1 are considered off)
					neighborStates.back().back().push_back(WireBrain->fullStateBuffer[neighborID].focalStates[j] == -1);
				}
				neighborStates.back().back().push_back(WireBrain->fullStateBuffer[i].focalStates[j] == 0); // add self, if not 0(WIRE STATE), since a node can't turn itself on, the only way it can effect future state is if it can or cannot be turned on!
                // so a strong self loop means that when this cell recives change, it's own state often inhibits the change (i.e., it's in charge or reset and thus ignores the neghbors state)
			}
			// now we can predict this focal states predictions!  We are using the neighbors (and self wire state) to predict the t+1 state so...
			std::cout << "running data flow for cell " << i << std::endl;
			auto trimmed_cellStates = TS::trimTimeSeries(cellStates.back(), TS::Position::FIRST, lifeTimes, 1);
			auto trimmed_neighborStates = TS::trimTimeSeries(neighborStates.back(), TS::Position::LAST, lifeTimes, 1);

			FRAG::saveFragMatrix(trimmed_cellStates, trimmed_neighborStates, "vis_output/cellPredictions/cellPrediction_" + std::to_string(i) + ".py", "feature", { "cell_" + std::to_string(i) }, -1, 0);
			FileManager::closeFile("vis_output/cellPredictions/cellPrediction_" + std::to_string(i) + ".py");

			std::vector<TS::intTimeSeries> decomposed_cellStates(12);
			std::vector <TS::intTimeSeries> decomposed_neighborStates(12);

			for (int t = 0; t < 12; t++) {
				std::cout << " creating states for cell " << i << " step " << t << std::endl;
				int samples = cellStates.back().size();
				for (int update = t; update < samples; update += 12) {
					decomposed_cellStates[t].push_back(cellStates.back()[update]);
					decomposed_neighborStates[t].push_back(neighborStates.back()[update]);
				}
			}

			// now that we have the decomposed states the trick is to use the neighbors from t to predict the cells at t+1
			for (int t = 0; t < 12; t++) {
				// start by predicting 0 which means use neighbors t-1 this means we need to trim the first of cell for each lifetime since that is the unpredicted initial state
				// and the last off neighbors to keep the lenghts the same
				if (t == 0) {
					auto trimmed_decomposed_cellStates = TS::trimTimeSeries(decomposed_cellStates[0], TS::Position::FIRST, decomposedLifetimes, 1); // trim the first, this is the predictee
					auto trimmed_decomposed_neighborStates = TS::trimTimeSeries(decomposed_neighborStates[11], TS::Position::LAST, decomposedLifetimes, 1); // trim the last, this is the predictor

					FRAG::saveFragMatrix(trimmed_decomposed_cellStates, trimmed_decomposed_neighborStates, "vis_output/cellPredictions/cellPrediction_" + std::to_string(i) + "_step_" + std::to_string(t) + ".py", "feature", { "cell_" + std::to_string(i) }, -1, 0);
					FileManager::closeFile("vis_output/cellPredictions/cellPrediction_" + std::to_string(i) + "_step_" + std::to_string(t) + ".py");

				}
				else {
					FRAG::saveFragMatrix(decomposed_cellStates[t], decomposed_neighborStates[t-1], "vis_output/cellPredictions/cellPrediction_" + std::to_string(i) + "_step_" + std::to_string(t) + ".py", "feature", { "cell_" + std::to_string(i) }, -1, 0);
					FileManager::closeFile("vis_output/cellPredictions/cellPrediction_" + std::to_string(i) + "_step_" + std::to_string(t) + ".py");
				}
			}

		}
		std::cout << "done with cell data flow" << std::endl;












		std::cout << "inflatedWorldStates.size() = " << inflatedWorldStates.size() << std::endl;

		// now we try to predict each world state from brain state
		// try to generate a frag matrix... 
 
		// first we need to create states for brain locations (we don't need inputs since the cells linked to inputs are set)
		// we need to itterate over all cells for every time step and create a state for each time step
		TS::intTimeSeries brainStates(WireBrain->brainStateLists.size());

		std::string cellToPredictorMap = "";

		for (int i : WireBrain->wireCellsActive) {
			//std::cout << i << std::endl;
			for (int j = 0; j < WireBrain->fullStateBuffer[i].focalStates.size();j++) {
				//std::cout << "   " << j << std::endl;
				brainStates[j].push_back(WireBrain->fullStateBuffer[i].focalStates[j]==-1);
			}
			cellToPredictorMap += std::to_string(WireBrain->fullStateBuffer[i].ID) + ",";
		}
		cellToPredictorMap.pop_back();
		FileManager::openAndWriteToFile("viz_output/cellToPredictorMap.csv", "map,\"" + cellToPredictorMap + "\"", "name, values");

		std::cout << "processed brain states size: " << brainStates.size() << std::endl;
		// now, do we have lifetimes?!
		std::cout << "WireBrain->lifeTimes count: " << WireBrain->lifeTimes.size() << std::endl;
		for (auto lt : WireBrain->lifeTimes) {
			std::cout << lt << "  ";
		}
		std::cout << std::endl;

		//auto trimmedBrainStates = TS::trimTimeSeries(brainStates, TS::Position::FIRST, WireBrain->lifeTimes, 4 * 12); // 4 world updates * 12 brain updates
		auto trimmedBrainStates = TS::trimTimeSeries(brainStates, TS::Position::FIRST, WireBrain->lifeTimes, currentLargestN * 12); // 5 world updates * 12 brain updates



		std::vector<std::string> featureNames({ "t-0","t-1","t-2","t-3","t-4"});

		//////////////////////////////////////////////////////////////////////////
		// full frag!
		//////////////////////////////////////////////////////////////////////////
		if (1) {
			std::cout << "Generating rawR full frag matrix" << std::endl;
			FRAG::saveFragMatrix(inflatedWorldStates, trimmedBrainStates,  "vis_output/rawR_FragmentationMatrix.py", "feature", featureNames, 2, 0);
		}
		//////////////////
		// decomposed frag
		//////////////////
		if (1) {

			int n = trimmedBrainStates.size() / 12;
			std::vector<std::vector<std::vector<int>>> brainStatesDecomposed(12);  // Create a vector with 12 groups
			std::vector<std::vector<std::vector<int>>> worldStatesDecomposed(12);  // Create a vector with 12 groups

			for (int i = 0; i < 12; ++i) {
				brainStatesDecomposed[i].resize(n);  // Each of the 12 groups will have n elements
				worldStatesDecomposed[i].resize(n);  // Each of the 12 groups will have n elements
				for (int j = 0; j < n; ++j) {
					brainStatesDecomposed[i][j] = trimmedBrainStates[i + 12 * j];  // Group every nth vector with offset by i
					worldStatesDecomposed[i][j] = inflatedWorldStates[i + 12 * j];  // Group every nth vector with offset by i
				}
			}

			if (0) { // switch on to check that decompose works
				for (int i = 0; i < 12 * 5; i++) {
					std::cout << i << ": ";
					for (const auto& val : inflatedWorldStates[i]) {
						std::cout << " " << val;
					}
					std::cout << std::endl;
				}

				for (int i = 0; i < 12; ++i) {
					std::cout << "Group " << i + 1 << ":\n";
					for (int j = 0; j < 5; ++j) {
						std::cout << "  Subgroup " << j + 1 << ":";
						for (const auto& val : worldStatesDecomposed[i][j]) {
							std::cout << " " << val;
						}
						std::cout << std::endl;
					}
				}

				for (int i = 0; i < 12 * 5; i++) {
					std::cout << i << ": ";
					for (const auto& val : trimmedBrainStates[i]) {
						std::cout << " " << val;
					}
					std::cout << std::endl;
				}

				for (int i = 0; i < 12; ++i) {
					std::cout << "Group " << i + 1 << ":\n";
					for (int j = 0; j < 5; ++j) {
						std::cout << "  Subgroup " << j + 1 << ":";
						for (const auto& val : brainStatesDecomposed[i][j]) {
							std::cout << " " << val;
						}
						std::cout << std::endl;
					}
				}

			}
			std::cout << worldStatesDecomposed[0].size() << "  " << inflatedWorldStates.size() << "  " << inflatedWorldStates[0].size() << std::endl;
			std::cout << brainStatesDecomposed[0].size() << "  " << trimmedBrainStates.size() << "  " << trimmedBrainStates[0].size() << std::endl;
			std::cout << "Generating time step rawR full frag matrix" << std::endl;
			for (int i = 0; i < 12; ++i) {
				FRAG::saveFragMatrix(worldStatesDecomposed[i], brainStatesDecomposed[i], "vis_output/rawR_step_" + std::to_string(i) + "_FragmentationMatrix.py", "feature", { "direction" }, 2, 0);
			}
		}




















		// now generate state to state information

		std::cout << "Generating State 2 State" << std::endl;

		// cast the brain to a wire cube brain so we can get at the guts!
		//auto WireBrain = std::dynamic_pointer_cast<WireCubeBrain>(groups[groupNameSpacePL->get(PT)]->population[0]->brains[brainNameSpacePL->get(PT)]);
		std::unordered_map<std::string, std::string> symbolMap; // for each rawState, it's mapped Name
		std::unordered_map<std::string, int> stateTransCounts; // for each transition, how often does it occur?
		std::unordered_map<std::string, std::vector<int>> stateTransWhen;// for each transition, when does it occur?

		std::string tempList;
		int nextSymbolCode = 0;
		std::string lastSymbol;
		double stateTransCountsMax = 0;
		int prtNxt = 0;
		int when = 0;
		for (auto& state : WireBrain->brainStateLists) {
			FileManager::openAndWriteToFile("vis_output/allBrainStates.txt", state);
			//std::cout << when << "   " << WireBrain->inputStatesStrings.size() << std::endl;
			//if ((when + 1) % lifeTimes[0] != 0) { // if this is the last step before a brain reset we should not save to state2state
			if (1) { // if this is the last step before a brain reset we should not save to state2state
				std::string in_and_state = WireBrain->inputStatesStrings[when] + "_" + state;
				//std::cout << in_and_state << std::endl;

				if (symbolMap.find(in_and_state) == symbolMap.end()) {
					// get new symbol
					int n = nextSymbolCode++;
					std::string symbol;
					while (n > 0) {
						//std::cout << n << std::endl;
						int remainder = n % 26;
						symbol = char('A' + remainder) + symbol;
						n /= 26;
					}
					if (symbol.empty()) symbol = "AAA";
					while (symbol.size() < 3) {
						symbol = "A" + symbol;
					}
					// get new symbol
					symbolMap[in_and_state] = symbol;
				}

				if (0) {
					if (symbolMap[in_and_state] == "AAA") {
						std::cout << symbolMap[in_and_state] << " > ";
						prtNxt = 60;
					}

					if (prtNxt > 0) {
						std::cout << symbolMap[in_and_state] << " > ";
						prtNxt--;
						if (prtNxt == 0) {
							std::cout << std::endl;
						}
					}
				}

				if (when % lifeTimes[0] > 0) {
					std::string transition = lastSymbol + "->" + symbolMap[in_and_state];
					stateTransCounts[transition]++;
					if (stateTransCounts[transition] > stateTransCountsMax) {
						stateTransCountsMax = stateTransCounts[transition];
					}
					if (stateTransWhen.find(transition) == stateTransWhen.end()) {
						stateTransWhen[transition] = std::vector<int>();
					}
					stateTransWhen[transition].push_back(when % lifeTimes[0]); // %lifeTimes[0] assumes all lifetimes are the same length!

				}
				lastSymbol = symbolMap[in_and_state];
			}
			when += 1;

		}

		std::string s2s_string = "digraph G {\n";

		for (auto& t : stateTransCounts) {
			double hue = 0;
			for (auto& v : stateTransWhen[t.first]) {
				hue += v;
			}
			s2s_string += t.first + " [color=\"" +
				std::to_string(.75 - (hue / double(stateTransWhen[t.first].size()) / double(lifeTimes[0] * (1.0 / .75)))) +
				",1,.8\",penwidth=" + std::to_string(3.0 + double(stateTransCounts[t.first]) / (stateTransCountsMax / 10.0)) + "]\n";
		}
		s2s_string += "}\n";
		FileManager::openAndWriteToFile("vis_output/state_2_state.txt", s2s_string);
	}

	////////////////////////////////////////////////////////
	// END OF NEW CODE TO SUPPORT WIRE BRAIN MEMORY PROJECT
	////////////////////////////////////////////////////////


	if (0) {
		auto lifeTimes = brain->getLifeTimes();
		auto inputStates = TS::remapToIntTimeSeries(brain->getInputStates(), TS::RemapRules::TRIT);
		auto outputStates = TS::remapToIntTimeSeries(brain->getOutputStates(), TS::RemapRules::TRIT);
		auto brainStates = TS::remapToIntTimeSeries(brain->getHiddenStates(), TS::RemapRules::TRIT);

		TS::intTimeSeries shortInputStates = TS::trimTimeSeries(inputStates, TS::Position::FIRST, lifeTimes, currentLargestN);

		TS::intTimeSeries shortOutputStatesBefore;// only needed if recurrent
		TS::intTimeSeries shortOutputStatesAfter; // always needed
		if (brain->recurrentOutput) {
			shortOutputStatesBefore = TS::trimTimeSeries(outputStates, TS::Position::FIRST, lifeTimes, currentLargestN); // trim off first currentLargestN for ramp up time
			auto tempStorterLifeTimes = TS::updateLifeTimes(lifeTimes, -1 * currentLargestN); // remove currentLargestN from each lifetime
			shortOutputStatesBefore = TS::trimTimeSeries(shortOutputStatesBefore, TS::Position::LAST, tempStorterLifeTimes, 1); // remove 1 from end
			shortOutputStatesAfter = TS::trimTimeSeries(outputStates, TS::Position::FIRST, lifeTimes, currentLargestN + 1); // remove currentLargestN+1 from front for ramp up + 1
		}
		else {
			shortOutputStatesAfter = TS::trimTimeSeries(outputStates, TS::Position::FIRST, lifeTimes, currentLargestN);
		}

		// always recurrent
		TS::intTimeSeries shortBrainStatesBefore = TS::trimTimeSeries(brainStates, TS::Position::FIRST, lifeTimes, currentLargestN);
		auto tempStorterLifeTimes = TS::updateLifeTimes(lifeTimes, -1 * currentLargestN);
		shortBrainStatesBefore = TS::trimTimeSeries(shortBrainStatesBefore, TS::Position::LAST, tempStorterLifeTimes, 1);
		TS::intTimeSeries shortBrainStatesAfter = TS::trimTimeSeries(brainStates, TS::Position::FIRST, lifeTimes, currentLargestN + 1);

		std::vector<int> shortLifeTimes = TS::updateLifeTimes(lifeTimes, -1 * currentLargestN);

		double R = ENT::ConditionalMutualEntropy(worldStates, shortBrainStatesAfter, shortInputStates);
		org->dataMap.append("R", R * RMultPL->get(PT));

		double rawR = ENT::MutualEntropy(worldStates, shortBrainStatesAfter);
		org->dataMap.append("rawR", rawR);





		double earlyRawR50 = ENT::MutualEntropy(TS::trimTimeSeries(worldStates, { 0,.5 }, shortLifeTimes), TS::trimTimeSeries(shortBrainStatesAfter, { 0,.5 }, shortLifeTimes));
		org->dataMap.append("earlyRawR50", earlyRawR50);

		double earlyRawR20 = ENT::MutualEntropy(TS::trimTimeSeries(worldStates, { 0,.2 }, shortLifeTimes), TS::trimTimeSeries(shortBrainStatesAfter, { 0,.2 }, shortLifeTimes));
		org->dataMap.append("earlyRawR20", earlyRawR20);

		double lateRawR50 = ENT::MutualEntropy(TS::trimTimeSeries(worldStates, { .5,1 }, shortLifeTimes), TS::trimTimeSeries(shortBrainStatesAfter, { .5,1 }, shortLifeTimes));
		org->dataMap.append("lateRawR50", lateRawR50);

		double lateRawR20 = ENT::MutualEntropy(TS::trimTimeSeries(worldStates, { .8,1 }, shortLifeTimes), TS::trimTimeSeries(shortBrainStatesAfter, { .8,1 }, shortLifeTimes));
		org->dataMap.append("lateRawR20", lateRawR20);






		if (analyze) {
			std::cout << "NBack World analyze... organism with ID " << org->ID << " scored " << org->dataMap.getAverage("score") << std::endl;
			FileManager::writeToFile("score_id_" + std::to_string(org->ID) + ".txt", std::to_string(org->dataMap.getAverage("score")));



			if (saveFragOverTime) { // change to 1 to save frag over time
				std::cout << "  saving frag over time..." << std::endl;
				std::string header = "LOD_order, score,";
				std::string outStr = std::to_string(org->dataMap.getIntVector("ID")[0]) + "," + std::to_string(org->dataMap.getAverage("score")) + ",";
				std::vector<int> save_levelsThresholds = { 50,75,100 };
				for (auto th : save_levelsThresholds) {
					auto frag = FRAG::getFragmentationSet(worldStates, shortBrainStatesAfter, ((double)th) / 100.0, "feature");
					for (int f = 0; f < frag.size(); f++) {
						header += "Threshold_" + std::to_string(th) + "__feature_" + std::to_string(f) + ",";
						outStr += std::to_string(frag[f]) + ",";
					}
				}
				FileManager::writeToFile("fragOverTime.csv", outStr.substr(0, outStr.size() - 1), header.substr(0, header.size() - 1));
			}

			if (saveBrainStructureAndConnectome) {
				std::cout << "saving brain connectome and structrue..." << std::endl;

				brain->saveConnectome("brainConnectome_id_" + std::to_string(org->ID) + ".py");
				brain->saveStructure("brainStructure_id_" + std::to_string(org->ID) + ".dot");
			}

			//if (1) {
			//	auto smearPair = SMR::getSmearednessConceptsNodesPair(shortInputStates, worldStates, shortBrainStatesAfter);
			//	FileManager::writeToFile("smear.txt", std::to_string(smearPair.second));
			//	FileManager::writeToFile("smear.txt", std::to_string(smearPair.first));
			//
			//	std::cout << "rawR for worldStateSet { i,j }, brainStateSet, { i,j } " << std::endl;
			//
			//	for (double i = 0; i <= 1; i += .1) {
			//		std::cout << i << " : ";
			//		for (double j = i + .1; j <= 1; j += .1) {
			//			std::cout << ENT::MutualEntropy(TS::trimTimeSeries(worldStates, { i,j }, shortLifeTimes), TS::trimTimeSeries(shortBrainStatesAfter, { i,j }, shortLifeTimes)) / ENT::Entropy(TS::trimTimeSeries(shortBrainStatesAfter, { i,j }, shortLifeTimes)) << " , ";
			//		}
			//		std::cout << std::endl;
			//	}
			//}


			if (saveStateToState) {
				std::cout << "  saving state to state..." << std::endl;
				std::string fileName = "StateToState_id_" + std::to_string(org->ID) + ".txt";
				if (brain->recurrentOutput) {
					S2S::saveStateToState({ brainStates, outputStates }, { inputStates }, lifeTimes, "H_O__I_" + fileName);
					S2S::saveStateToState({ brainStates }, { inputStates }, lifeTimes, "H_I_" + fileName);
				}
				else {
					S2S::saveStateToState({ brainStates, TS::extendTimeSeries(outputStates, lifeTimes, {0}, TS::Position::FIRST) }, { inputStates }, lifeTimes, "H_O__I_" + fileName);
					S2S::saveStateToState({ brainStates }, { outputStates, inputStates }, lifeTimes, "H__O_I_" + fileName);
					S2S::saveStateToState({ brainStates }, { inputStates }, lifeTimes, "H_I_" + fileName);
				}
			}
			//std::cout << "worldEnt: " << ENT::Entropy(worldStates) << "  brainEnt: " << ENT::Entropy(shortBrainStatesAfter) << "  worldBrainEnt: " << ENT::Entropy(TS::Join(worldStates, shortBrainStatesAfter)) << "  rawR: " << rawR << std::endl;
			//std::cout << "earlyRawR20: " << earlyRawR20 << "  earlyRawR50: " << earlyRawR50 << "  lateRawR50: " << lateRawR50 << "  lateRawR20: " << lateRawR20 << std::endl;

			// save fragmentation matrix of brain(hidden) predictions of world features
			if (save_R_FragMatrix) {
				std::cout << "  saving R frag matrix..." << std::endl;

				FRAG::saveFragMatrix(worldStates, shortBrainStatesAfter, "R_FragmentationMatrix_id_" + std::to_string(org->ID) + ".py", "feature");
			}
			// save data flow information - 
			//std::vector<std::pair<double, double>> flowRanges = { {0,1},{0,.333},{.333,.666},{.666,1},{0,.5},{.5,1} };
			////std::vector<std::pair<double, double>> flowRanges = { {0,1},{.5,1} };///, { 0,.1 }, { .9,1 }};


			if (saveFlowMatrix) {
				std::cout << "  saving flow matix..." << std::endl;
				std::vector<std::pair<double, double>> flowRanges = { {0,.25}, {.75,1}, {0,1} };//, { 0,.1 }, { .9,1 }};

				//std::cout << TS::TimeSeriesToString(TS::trimTimeSeries(brainStates, TS::Position::LAST, lifeTimes), ",",",") << std::endl;
				//std::cout << TS::TimeSeriesToString(TS::trimTimeSeries(brainStates, TS::Position::FIRST, lifeTimes), ",",",") << std::endl;
				if (brain->recurrentOutput) {
					FRAG::saveFragMatrixSet(
						TS::Join({ TS::trimTimeSeries(brainStates, TS::Position::FIRST, lifeTimes), TS::trimTimeSeries(outputStates, TS::Position::FIRST, lifeTimes) }),
						TS::Join({ TS::trimTimeSeries(brainStates, TS::Position::LAST, lifeTimes), inputStates, TS::trimTimeSeries(outputStates, TS::Position::LAST, lifeTimes) }),
						lifeTimes, flowRanges, "flowMap_id_" + std::to_string(org->ID) + ".py", "shared", -1);
				}
				else {
					FRAG::saveFragMatrixSet(
						TS::Join(TS::trimTimeSeries(brainStates, TS::Position::FIRST, lifeTimes), outputStates),
						TS::Join(TS::trimTimeSeries(brainStates, TS::Position::LAST, lifeTimes), inputStates),
						lifeTimes, flowRanges, "flowMap_id_" + std::to_string(org->ID) + ".py", "shared", -1);
				}
			}
			//auto flowMatrix = FRAG::getFragmentationMatrix(TS::Join(TS::trimTimeSeries(brainStates, TS::Position::FIRST, lifeTimes), outputStates), TS::Join(TS::trimTimeSeries(brainStates, TS::Position::LAST, lifeTimes), inputStates), "feature");




			if (saveStates) {
				std::cout << "  saving brain states information..." << std::endl;
				std::string fileStr = "";
				if (brain->recurrentOutput) {

					auto discreetInput = inputStates;
					auto discreetOutputBefore = TS::trimTimeSeries(outputStates, TS::Position::LAST, lifeTimes);;
					auto discreetOutputAfter = TS::trimTimeSeries(outputStates, TS::Position::FIRST, lifeTimes);
					auto discreetHiddenBefore = TS::trimTimeSeries(brainStates, TS::Position::LAST, lifeTimes);
					auto discreetHiddenAfter = TS::trimTimeSeries(brainStates, TS::Position::FIRST, lifeTimes);

					fileStr += "input,outputBefore,outputAfter,hiddenBefore,hiddenAfter\n";
					for (int i = 0; i < discreetInput.size(); i++) {
						fileStr += "\"" + TS::TimeSeriesSampleToString(discreetInput[i], ",") + "\",";
						fileStr += "\"" + TS::TimeSeriesSampleToString(discreetOutputBefore[i], ",") + "\","; // every other
						fileStr += "\"" + TS::TimeSeriesSampleToString(discreetOutputAfter[i], ",") + "\","; // the other ones
						fileStr += "\"" + TS::TimeSeriesSampleToString(discreetHiddenBefore[i], ",") + "\","; // every other
						fileStr += "\"" + TS::TimeSeriesSampleToString(discreetHiddenAfter[i], ",") + "\"\n"; // the other ones
					}
				}
				else {

					auto discreetInput = inputStates;
					auto discreetOutput = outputStates;
					auto discreetHiddenBefore = TS::trimTimeSeries(brainStates, TS::Position::LAST, lifeTimes);
					auto discreetHiddenAfter = TS::trimTimeSeries(brainStates, TS::Position::FIRST, lifeTimes);

					fileStr += "input,output,hiddenBefore,hiddenAfter\n";
					for (int i = 0; i < discreetInput.size(); i++) {
						fileStr += "\"" + TS::TimeSeriesSampleToString(discreetInput[i], ",") + "\",";
						fileStr += "\"" + TS::TimeSeriesSampleToString(discreetOutput[i], ",") + "\",";
						fileStr += "\"" + TS::TimeSeriesSampleToString(discreetHiddenBefore[i], ",") + "\","; // every other
						fileStr += "\"" + TS::TimeSeriesSampleToString(discreetHiddenAfter[i], ",") + "\"\n"; // the other ones
					}
				}

				FileManager::writeToFile("NBack_BrainActivity_id_" + std::to_string(org->ID) + ".csv", fileStr);


			}
			std::cout << "  ... analyze done" << std::endl;
		} // end analyze

	}


}

std::unordered_map<std::string, std::unordered_set<std::string>> NBackWorld::requiredGroups() {
	return { {groupName, {"B:" + brainName + ",1," + std::to_string(N2OutMap.size())}} };
}



