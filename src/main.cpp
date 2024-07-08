#define HAVE_SSTREAM

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <vector>
#include <random>
#include <limits>
#include <algorithm>
#include <cmath>
#include <thread>
#include <mutex>

#include "CSVHandler.h"


// the general include for eo
#include <eo>
#include <es.h>
#include <eval/moFullEvalByModif.h>
#include <eval/moFullEvalByCopy.h>
#include <eoInt.h>
#include <utils/eoParam.h>

// Fitness function
#include <eoEvalFunc.h>

// Neighbors and Neighborhoods
#include <neighborhood/moNeighbor.h>
#include <neighborhood/moNeighborhood.h>

// ########## Algorithm and its components ##########
// Hill-climber and it components
#include <algo/moSimpleHC.h>
// Simulated Annealing
#include <algo/moSA.h>

// Cooling Schedule
#include <coolingSchedule/moCoolingSchedule.h>

// Comparator
#include <comparator/moSolNeighborComparator.h>

//continuators
#include <continuator/moTrueContinuator.h>
#include <continuator/moCheckpoint.h>
#include <continuator/moCounterStat.h>
#include <continuator/moTimeContinuator.h>

// Tabu search
#include <algo/moTS.h>


//-----------------------------------------------------------------------------
// Global problem description
// The number of agents, and the number of items
thread_local int N, M;

// Valuations[n][m] is the utility of agent n and item m
thread_local std::vector<std::vector<int>> Valuations;

//-----------------------------------------------------------------------------
// representation of solutions, and neighbors
// eoInt is based on eoVector class
// Change <unsigned int> for maximin and <double> with leximin
typedef eoInt<unsigned int> Allocation;

//#########################################################
// REPRESENTATION
// define Population of solutions
static std::vector<Allocation> Population;
//#########################################################

//-----------------------------------------------------------------------------
// fitness function and evaluation of solution
template< class EOT >
class allocationEval : public eoEvalFunc<EOT>
{
public:

    /**
     * Find the agent with min utility
     */            
    void operator()(EOT& allocation) {
        static thread_local std::vector<typename EOT::Fitness> utility;

        // Set each agent's evaluation to zero
        utility.resize(N);
        std::fill(utility.begin(), utility.end(), 0);

        // For each item, increase the utility of the agent who gets the item
        for (int m = 0; m < allocation.size(); m++) {
            int agent = allocation[m];
            utility[agent] += Valuations[agent][m];
        }

        // sort utilitites in ascending order. This is used for leximin calculations
        std::sort(utility.begin(), utility.end());

        /*
        // Debug output to see the utilities
        std::cout << "Utilities: ";
        for (auto u : utility) {
            std::cout << u << " ";
        }
        std::cout << std::endl; */

        // calculate fitnes with leximin
        
        // TODO: endre 0.1 til 1/(maxUti + 1)
        /*
        const double k = 0.01; 
        typename EOT::Fitness fitness = 0;
        
        
        // iterate from smallest to largest utility
        for (std::size_t i = 0; i < utility.size(); i++)
        {
            fitness += utility[i] * std::pow(k, i);
        }*/
        

        // Maximin allocation
        
        // Find the agent with the lowest utility
        auto minUtility = std::numeric_limits<typename EOT::Fitness>::max();
        for (auto util : utility)
            minUtility = std::min(minUtility, util);

        // Inform the allocation about its utility
        // fitness with leximin and minUtility with maximin
        allocation.fitness(minUtility);

        // Debug output to confirm fitness value
        //std::cout << "Computed fitness: " << fitness << std::endl;
    }
};

//-----------------------------------------------------------------------------
// neighbor description
// Neighbor = How to compute the neighbor from the solution + information on it (i.e. fitness)
// all classes from paradisEO-mo use this template type

/**
 * Modifying an existing allocation by moving a given item to a given agent
 */
//template<class EOT, class Fitness = unsigned int>
template<class EOT, class Fitness = typename EOT::Fitness>
class moveNeighbor: public moNeighbor<EOT, Fitness> {
public:

	/**
	 * Apply the move to a given solution
	 * @param _solution the solution to move
	 */
	virtual void move(EOT& solution) {
        //std::cout << "Applying move: " << item << " to agent " << agent << std::endl;
        solution[item] = agent;
        solution.invalidate();
	}

	void setMove(int item, int agent) {
		this->item = item;
        this->agent = agent;
	}

	void getMove(int & item, int & agent) {
        item = this->item;
        agent = this->agent;
	}

	virtual bool equals(moveNeighbor<EOT,Fitness>& neighbor) {
		return item == neighbor.item && agent == neighbor.agent;
	}

	void print() {
		std::cout << "[give " << item << " to " << agent << "] -> " << this->fitness() << std::endl;
	}

private:
    int item;
    int agent;
};

//-----------------------------------------------------------------------------
// neighborhood description
// template <class EOT, class Fitness= unsigned int>
template <class EOT, class Fitness=typename EOT::Fitness>
class moveNeighborhood : public moNeighborhood<moveNeighbor<EOT, Fitness> >
{
public:
    typedef moveNeighbor<EOT, Fitness> Neighbor;

    std::vector<std::pair<int, int>> allMoves;
    typename std::vector<std::pair<int, int>>::iterator currentMove;

    moveNeighborhood() {
        for (int item = 0; item < M; ++item) {
            for (int agent = 0; agent < N; ++agent) {
                allMoves.emplace_back(item, agent);
            }
        }
    }

    /**
     * @return true if there is at least an available neighbor
     */
    virtual bool hasNeighbor(EOT& solution) {
        return !allMoves.empty();
    }

    /**
     * @return if neighborhood is random (default false)
    */
    virtual bool isRandom() override {
        return true;
    }

    /**
     * Initialization of the neighborhood
     * @param solution the solution to explore
     * @param current the first neighbor
     */
    virtual void init(EOT& solution, Neighbor& current) {
        std::shuffle(allMoves.begin(), allMoves.end(), std::default_random_engine());
        currentMove = allMoves.begin();
        current.setMove(currentMove->first, currentMove->second);
    }

    /**
     * Give the next neighbor
     * @param solution the solution to explore
     * @param current the next neighbor
     */
    virtual void next(EOT& solution, Neighbor& current) {
        ++currentMove;
        if (currentMove != allMoves.end()) {
            current.setMove(currentMove->first, currentMove->second);
        }
    }

    
    virtual bool cont(EOT& _solution) {
        return currentMove != allMoves.end();
    }

    virtual std::string className() const {
        return "moveNeighborhood";
    }
};

// ##################################################################################
// ##################   HELP FUNCTIONS FOR GENETIC ALGORITHMS #######################

// Function to generate a single solution
Allocation generateSolution( int m, std::mt19937& gen) {
    std::uniform_int_distribution<> distr(0, N-1);
    Allocation solution(N);
    for (int i = 0; i < N; i++) {
        solution[i] = distr(gen);
    }
    return solution;
}

// Function to generate multiple solutions
std::vector<Allocation> generateMultipleSolutions(int numSolutions, int m) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<Allocation> solutions;
    solutions.reserve(numSolutions);

    for (int i = 0; i < numSolutions; i++) {
        solutions.push_back(generateSolution(m, gen));
    }
    return solutions;
}

// Function to evaluate and print solutions
void evaluateAndPrintSolutions(std::vector<Allocation>& solutions) {
    allocationEval<Allocation> fullEval; 

    for (size_t i = 0; i < solutions.size(); i++) {
        fullEval(solutions[i]);  // Evaluate the solution
        
        // Print the solution
        std::cout << "Solution " << i+1 << ":" << std::endl;
        std::cout << solutions[i] << std::endl << std::endl;
    }
}

// Helper function
std::vector<int> convertToVectorOfInt(const Allocation& allocation) {
    std::vector<int> intVector;
    for (const auto& value : allocation) {
        intVector.push_back(static_cast<int>(std::round(value))); // Convert and round the double to int
    }
    return intVector;
}

void run_algorithms (int thread_local_N, const std::string& filepath_executions, const std::string& filepath_valuations, const std::string& filepath_solutions) 
{
    // Set thread-local N to the passed value
    N = thread_local_N;

    // Initialization of results vector
    std::vector<std::pair<std::string, std::vector<int>>> results = {
        {"HC", {}},
        {"HC01", {}},
        {"SA100", {}},
        {"SA5", {}},
        {"Number of items", {}},
        {"Profile ID", {}}
    };

    // To store all valuation profiles from iterations
    std::vector<std::vector<std::vector<int>>> allValuations;

    std::vector<std::vector<std::vector<int>>> finalSolution;
    

    /* =========================================================
     *
     * the cooling schedule of the process
     *
     * ========================================================= */
    // Define the parameters for the cooling schedule
    double initialTemperature = 1;      // Initial temperature
    double coolingFactor = 0.9;        // Factor by which the temperature will be multiplied at each step
    unsigned span = 100;              // Number of iterations with the same temperature
    unsigned s = 5;                 // Number of iterations with the same temperature
    double finalTemperature = 0.01; // Final temperature, stopping criterion

    moSimpleCoolingSchedule<Allocation> coolingSchedule(initialTemperature, coolingFactor, span, finalTemperature);
    moSimpleCoolingSchedule<Allocation> divtwospan(initialTemperature, coolingFactor, span, finalTemperature);

    /* =========================================================
     *
     * Checkpointing
     *
     * ========================================================= */
    moTrueContinuator<moveNeighbor<Allocation>> continuator;
    moCheckpoint<moveNeighbor<Allocation>> checkpoint(continuator);
    moCounterStat<Allocation> iterStat;
    checkpoint.add(iterStat);
    moTimeContinuator<moveNeighbor<Allocation>> milli(0.001, false);
    moTimeContinuator<moveNeighbor<Allocation>> smaller(0.0001, false);

    /* =========================================================
     *
     * Initialization and Execution of algorithms
     *
     * ========================================================= */

    for (int i = 0; i < 100; i++)
    {   
        // Create a random device and seed a generator
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distr(0,N-1);
        // Generate a random number between 2*N and 4*N
        std::uniform_int_distribution<int> uni(2 * N, 4 * N);
        M = uni(gen);
        std::cout << "Number of items (M): " << M << std::endl;

        // Generate single valuation vector for items
        /*
        std::vector<int> single_valuation(M);
        for (int m = 0; m < M; m++)
        {
            single_valuation[m] = distr(gen);
        }
        
        // Create Valuation matrix where each agent has the same valuations
        Valuations.clear();
        Valuations.resize(N, single_valuation);*/
        /*
        // Print the Valuations matrix
        for (int n = 0; n < N; ++n) {
            for (int m = 0; m < M; ++m) {
                std::cout << Valuations[n][m] << " ";
            }
            std::cout << std::endl;
        }*/

        
        // Random valuation profile
        
        Valuations.clear();
        Valuations.resize(N, std::vector<int>(M));
        // Assign random valuations to each item
        for (int n = 0; n < N; n++) {
            for (int m = 0; m < M; m++) {
                Valuations[n][m] = distr(gen);
            }
        }

        // Store current Valuations matrix to allValuations
        allValuations.push_back(Valuations);

        // Printing out the Valuations matrix
        /*
        std::cout << "Valuations Matrix:" << std::endl;
        for (int n = 0; n < N; n++) {
            for (int m = 0; m < M; m++) {
                std::cout << Valuations[n][m] << " ";
            }
            std::cout << std::endl; // New line at the end of each row
        }*/
        
        // Initialize solution with random elements from 1 to n
        Allocation solution(M);
        for (int m = 0; m < M; m++)
        {
            solution[m] = distr(gen);
        }
        
        Allocation tmp1 = solution;
        Allocation tmp2 = solution;
        Allocation tmp3 = solution; 
        Allocation tmp4 = solution; 

        // EVAL INSTANCES
        // Create eval instances
        allocationEval<Allocation> fullEvalHc;
        allocationEval<Allocation> fullEvalSa;
        allocationEval<Allocation> fullEvalTs;
        moFullEvalByCopy<moveNeighbor<Allocation>> moveEvalHC(fullEvalHc);
        moFullEvalByCopy<moveNeighbor<Allocation>> moveEvalSa(fullEvalSa);
        moFullEvalByCopy<moveNeighbor<Allocation>> moveEvalTs(fullEvalTs);

        // Evaluatig initial solution
        //fullEvalHc(tmp1);
        //fullEvalHc(tmp2);
        //fullEvalSa(tmp3);
        //fullEvalSa(tmp4);
        fullEvalTs(tmp1);
        fullEvalTs(tmp2);
        fullEvalTs(tmp3);
        fullEvalTs(tmp4);

        //move neighborhood
        moveNeighborhood<Allocation> moveNH;

        
        //std::cout << "Initial Solution:" << std::endl;
        //std::cout << solution << std::endl << std::endl;
        //std::cout << tmp1 << std::endl << std::endl;
        //std::cout << tmp2 << std::endl << std::endl;
    

        // EXECUTING algorithms
        //moSimpleHC<moveNeighbor<Allocation>> hc(moveNH, fullEvalHc, moveEvalHC, milli); // removed last paramater milli, but that backlater
        //hc(tmp1);
        //moSimpleHC<moveNeighbor<Allocation>> hc_smaller(moveNH, fullEvalHc, moveEvalHC, smaller);
        //hc_smaller(tmp2);

        //moSA<moveNeighbor<Allocation>> sa(moveNH, fullEvalSa, moveEvalSa, coolingSchedule, milli);
        //sa(tmp3);
        //moSA<moveNeighbor<Allocation>> sa_s(moveNH, fullEvalSa, moveEvalSa, divtwospan, milli);
        //sa_s(tmp4);

        
        moTS<moveNeighbor<Allocation>> ts_sec(moveNH, fullEvalTs, moveEvalTs, 1, 20);
        moTS<moveNeighbor<Allocation>> ts_ms(moveNH, fullEvalTs, moveEvalTs, 0.01, 20);
        moTS<moveNeighbor<Allocation>> ts_big(moveNH, fullEvalTs, moveEvalTs, 0.001, 20);
        moTS<moveNeighbor<Allocation>> ts_smal(moveNH, fullEvalTs, moveEvalTs, 0.0001, 20);
        ts_sec(tmp1);
        ts_ms(tmp2);
        ts_big(tmp3);
        ts_smal(tmp4);

        // Create a 2D vector and add the Allocation object
        std::vector<int> vec1D = convertToVectorOfInt(tmp1);

        // Create a 2D vector and add the 1D vector
        std::vector<std::vector<int>> vec2D;
        vec2D.push_back(vec1D);

        // Add the 2D vector to the 3D vector
        finalSolution.push_back(vec2D);

        
        //std::cout << "tmp1: " << tmp1 << std::endl << std::endl ;
        //std::cout << "tmp2: " << tmp2 << std::endl << std::endl ;
        //std::cout << "tmp3: " << tmp3 << std::endl << std::endl ;
        //std::cout << "tmp4: " << tmp4 << std::endl << std::endl ;
        //std::cout << "final TS: " << tmp3 << std::endl << std::endl ;
        

        // Add fitness values to vectors in results
        results[0].second.push_back(tmp1.fitness());
        results[1].second.push_back(tmp2.fitness());
        results[2].second.push_back(tmp3.fitness());
        results[3].second.push_back(tmp4.fitness());
        results[4].second.push_back(1.0 * M);
        results[5].second.push_back(1.0 * i);
    }

    // Writing stored results and their valuation-profiles to CSV
    {   
        CSVHandler::write_matrix_csv(filepath_solutions, finalSolution);
        CSVHandler::write_matrix_csv(filepath_valuations, allValuations);
        CSVHandler::write_csv(filepath_executions, results);
    }

}


void gen_alg (int thread_local_N, const std::string& filepath_executions, const std::string& filepath_valuations)
{
    // Set thread-local N to the passed value
    N = thread_local_N;

    // Initialization of results vector
    std::vector<std::pair<std::string, std::vector<int>>> results = {
        {"GA Solution", {}}
    };

    // To store all valuation profiles from iterations
    std::vector<std::vector<std::vector<int>>> allValuations;

    for (int i = 0; i < 100; i++)
    {
        // Create a random device and seed a generator
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distr(0,N-1);
        // Generate a random number between 2*N and 4*N
        std::uniform_int_distribution<int> uni(2 * N, 4 * N);
        M = uni(gen);
        //std::cout << "Number of items (M): " << M << std::endl;

        Valuations.clear();
        Valuations.resize(N, std::vector<int>(M));
        // Assign random valuations to each item
        for (int n = 0; n < N; n++) {
            for (int m = 0; m < M; m++) {
                Valuations[n][m] = distr(gen);
            }
        }

        // Store current Valuations matrix to allValuations
        allValuations.push_back(Valuations);

        // Printing out the Valuations matrix
        /*
        std::cout << "Valuations Matrix:" << std::endl;
        for (int n = 0; n < N; n++) {
            for (int m = 0; m < M; m++) {
                std::cout << Valuations[n][m] << " ";
            }
            std::cout << std::endl; // New line at the end of each row
        }*/

        // Create eval instances
        allocationEval<Allocation> fullEval;
        moFullEvalByCopy<moveNeighbor<Allocation>> moveEval(fullEval);
        
        const unsigned int T_SIZE = 8; // size for tournement selection
        const unsigned int VEC_SIZE = M; // Number of object variables in genotypes (number of elements in the vector)
        const unsigned int POP_SIZE = 250; // Size of population

        const unsigned int MAX_GEN = 500000; // Maximum number of generation before STOP
        const unsigned int MIN_GEN = 10000;  // Minimum number of generation before STOP
        const unsigned int STEADY_GEN = 10000;  //stop after STEADY_GEN gen. without improvement

        const float P_CROSS = 0.67;      // Crossover probability
        const float P_MUT = 0.1;        // Mutation probability 

        const double EPSILON = 0.01;	// range for real uniform mutation
        double SIGMA = 0.3;	    	// std dev. for normal mutation
        
        const double hypercubeRate = 0.6;     // relative weight for hypercube Xover
        const double segmentRate = 0.6;  // relative weight for segment Xover

        const double uniformMutRate = 0.3;  // relative weight for uniform mutation
        const double detMutRate = 0.3;      // relative weight for det-uniform mutation
        const double normalMutRate = 0.3;   // relative weight for normal mutation

        // ############################
        // Initialisation of population
        // ############################
        eoUniformGenerator<size_t> uGen(0, N); // hadde N-1 her
        eoInitFixedLength<Allocation> random(VEC_SIZE, uGen);
        // Initialization of the population
        eoPop<Allocation> pop(POP_SIZE, random);

        // and evaluate it in one loop
        apply<Allocation>(fullEval, pop);

        // OUTPUT
        pop.sort();
        
        std::cout << "Initial Population" << std::endl;
        std::cout << pop.best_element() << std::endl;

        // ############################
        //          Selection
        // ############################
        //eoDetTournamentSelect<Allocation> selectOne(T_SIZE);
        eoProportionalSelect<Allocation> selectOne(pop);
        eoSelectPerc<Allocation> select(selectOne);

        eoGenerationalReplacement<Allocation> replace;
        eoWeakElitistReplacement<Allocation> addWeakElitisim(replace); // la inn addWeak.. istedetfor kun "replace" lenger ned

        //eoPlusReplacement<Allocation> replace;
        eoSegmentCrossover<Allocation> xoverS;

        eoHypercubeCrossover<Allocation> xoverA;

        eoPropCombinedQuadOp<Allocation> xover(xoverS, segmentRate);
        xover.add(xoverA, hypercubeRate);

        // ############################
        //          MUTATION
        // ############################
        // offspring(i) uniformly chosen in [parent(i)-epsilon, parent(i)+epsilon]
        eoUniformMutation<Allocation>  mutationU(EPSILON);
        
        // k (=1) coordinates of parents are uniformly modified
        eoDetUniformMutation<Allocation>  mutationD(EPSILON);
        // all coordinates of parents are normally modified (stDev SIGMA)
        //eoNormalMutation<Allocation>  mutationN(SIGMA);

        
        // Combine them with relative weights
        eoPropCombinedMonOp<Allocation> mutation(mutationU, uniformMutRate);
        mutation.add(mutationD, detMutRate);
        //mutation.add(mutationN, normalMutRate, true);

        // ############################
        //  CHECKPOINT
        // ############################
        // termination conditions: use more than one

        //STOP after MAX_GEN generations
        eoGenContinue<Allocation> genCont(MAX_GEN);
        // do MIN_GEN gen., then stop after STEADY_GEN gen. without improvement
        eoSteadyFitContinue<Allocation> steadyCont(MIN_GEN, STEADY_GEN);
        // stop when fitness reaches a target (here VEC_SIZE)
        //eoSteadyFitContinue<Allocation> fitCont(0);

        // do stop when one of the above says so
        eoCombinedContinue<Allocation> gaContinuator(genCont);
        gaContinuator.add(steadyCont);

        // The operators are  encapsulated into an eoTRansform object
        eoSGATransform<Allocation> transform(xover, P_CROSS, mutation, P_MUT);

        // ############################
        //          GENERATION
        // ############################
        //         the algorithm

        eoEasyEA<Allocation> gga(gaContinuator, fullEval, select, transform, addWeakElitisim);
        // Apply algo on the population
        gga(pop);

        //OUTPUT
        
        //std::cout << "########################" << std::endl << std::endl;
        pop.sort();
        //std::cout << "FINAL Population\n" << pop << std::endl;

        //std::cout << "########################" << std::endl << std::endl;
        //std::cout << "The best final element: \n" << pop.best_element() << std::endl; 
        //std::cout << "########################" << std::endl << std::endl;
        //std::cout << "########################" << std::endl << std::endl;
        //std::cout << "########################" << std::endl << std::endl;
        //std::cout << "Next iteration starts now" << std::endl << std::endl;
        //std::cout << "########################" << std::endl << std::endl;
        //std::cout << "########################" << std::endl << std::endl;
        //std::cout << "########################" << std::endl << std::endl;

        results[0].second.push_back(pop.nth_element_fitness(0));
    }
    

    {
        CSVHandler::write_matrix_csv(filepath_valuations, allValuations);
        CSVHandler::write_csv(filepath_executions, results);
    }


}

//-----------------------------------------------------------------------------

int main()
{
    /*
    N = 10;
   
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0,N-1);
    std::uniform_int_distribution<int> uni(2 * N, 4 * N);
    // Generate a random number between 2*N and 4*N
    M = uni(gen);
    //std::cout << "Number of items (M): " << M << std::endl;

    // Generate single valuation vector for items
    std::vector<int> single_valuation(M);
    for (int m = 0; m < M; m++)
    {
        single_valuation[m] = distr(gen);
    }
        
    // Create Valuation matrix where each agent has the same valuations
    Valuations.clear();
    Valuations.resize(N, single_valuation);*/

    // Print the Valuations matrix
    /*
    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < M; ++m) {
            std::cout << Valuations[n][m] << " ";
        }
        std::cout << std::endl;
    }*/

    /*
    Valuations.clear();
    Valuations.resize(N, std::vector<int>(M));
    // Assign random valuations to each item
    
    for (int n = 0; n < N; n++) {
        for (int m = 0; m < M; m++) {
            Valuations[n][m] = distr(gen);
        }
    }

    // Printing out the Valuations matrix
        
    std::cout << "Valuations Matrix:" << std::endl;
    for (int n = 0; n < N; n++) {
        for (int m = 0; m < M; m++) {
            std::cout << Valuations[n][m] << " ";
        }
        std::cout << std::endl; // New line at the end of each row
    }*/

    // ####### Manually setup for testing purposes #######
    // Manually set valuations for the agents

    //Valuations[0] = {1, 2, 3, 4, 5, 6, 7, 4}; // First agent
    //Valuations[1] = {4, 3, 2, 1, 2, 5, 6, 7}; // Second agent
    //Valuations[2] = {2, 1, 1, 4, 1, 1, 7, 7}; // Third agent
    //Valuations[3] = {7, 7, 7, 7, 7, 0, 0, 1}; // Fourth agent
    //Valuations[4] = {0, 0, 0, 6, 6, 6, 7, 7}; // Fifth agent
    //Valuations[5] = {7, 6, 5, 4, 3, 2, 1, 0}; // Sixth agent
    //Valuations[6] = {1, 1, 3, 1, 5, 6, 6, 7}; // Seventh agent
    
   
    /*
    // Initialize solution with random elements from 1 to n
    Allocation solution(M);
    for (int m = 0; m < M; m++)
    {
        solution[m] = distr(gen);
    }

    //Allocation tmp1 = solution; // For 10 seconds
    //Allocation tmp2 = solution; // For 1 ms

    // Create eval instances
    allocationEval<Allocation> fullEval;
    moFullEvalByCopy<moveNeighbor<Allocation>> moveEval(fullEval);

    // Evaluatig initial solution
    //fullEval(tmp1);
    //fullEval(tmp2);
    fullEval(solution);

    std::cout << "Initial Solution:" << std::endl;
    std::cout << solution << std::endl << std::endl;
    
    //move neighborhood
    moveNeighborhood<Allocation> moveNH;*/

    // Iterate through all neighbours and print their evaluation
    
    //int moveCounter = 0;
    

    /*
    moveNeighbor<Allocation> move;
    moveNH.init(solution, move);
    moveEval(solution, move);
    move.print();
    while(moveNH.cont(solution)) {
        moveNH.next(solution, move);
        moveEval(solution, move);
        move.print();
        moveCounter++;
    }

    std::cout << "moves: "<< moveCounter << std::endl <<std::endl;
    */

    /* =========================================================
     *
     *          the cooling schedule of the process
     *
     * ========================================================= */
    // Define the parameters for the cooling schedule
    /*
    double initialTemperature = 1;      // Initial temperature
    double coolingFactor = 0.9;        // Factor by which the temperature will be multiplied at each step
    unsigned span = 100;              // Number of iterations with the same temperature
    double finalTemperature = 0.01;  // Final temperature, stopping criterion

    moSimpleCoolingSchedule<Allocation> coolingSchedule(initialTemperature, coolingFactor, span, finalTemperature);*/

    /* =========================================================
     *
     * Checkpointing
     *
     * ========================================================= */
    /*
    moTrueContinuator<moveNeighbor<Allocation>> continuator;
    moCheckpoint<moveNeighbor<Allocation>> checkpoint(continuator);
    moCounterStat<Allocation> iterStat;
    checkpoint.add(iterStat);
    moTimeContinuator<moveNeighbor<Allocation>> t(0.01, false);
    checkpoint.add(t);*/


    /* =========================================================
     *
     * Initialization and Execution of algorithms
     *
     * ========================================================= */
     //Define the simple Hill-Climbing 
    /*
    moSimpleHC<moveNeighbor<Allocation>> hc(moveNH, fullEval, moveEval,t);
    // apply the local search on the solution
    hc(solution);

    // Simulated Annealing
    //moSA<moveNeighbor<Allocation>> localSearch1(moveNH, fullEval, moveEval, coolingSchedule, checkpoint);
    //localSearch1(solution);

    // Printing result of Simulated annealing
    std::cout << "final: " << solution << std::endl << std::endl ;
    std::cout << "Iterations: " << iterStat.getValue() << std::endl ;*/



    //Tabu search
    /*
    moTS<moveNeighbor<Allocation>> ls_sec(moveNH, fullEval, moveEval, 10, 10);
    moTS<moveNeighbor<Allocation>> ls_ms(moveNH, fullEval, moveEval, 2, 10);

    ls_sec(tmp1);
    ls_ms(tmp2);

    //Printing out results
    std::cout << "final: " << tmp1 << std::endl << std::endl ;
    //std::cout << "Iterations: " << iterStat.getValue() << std::endl ;
    std::cout << "#####################################" << std::endl;
    std::cout << "final: " << tmp2 << std::endl << std::endl ;*/


    /**=========================================================
     * 
     * 
     * EXECUTING MAIN PROGRAM
     * 
     *==========================================================*/
    
    // EXECUTE THREADS
    std::vector<std::thread> threads;
    std::vector<int> N_values = {2, 3, 4, 5, 6, 7, 8, 9, 10};
    //std::vector<int> N_values = {8};
    /*
    for (int N_value : N_values) {
        std::string filepath_executions = "lexgacsv/" + std::to_string(N_value) + "agents_executions.csv";
        std::string filepath_valuations = "lexgacsv/"+ std::to_string(N_value) + "agents_valuations.csv";
        threads.emplace_back(gen_alg, N_value, filepath_executions, filepath_valuations);
    }

    for (auto& thread : threads) {
        thread.join();
    }*/

    
    for (int N_value : N_values) {
        //std::string filepath_executions = "mmcsv/" + std::to_string(N_value) + "agents_executions.csv";
        //std::string filepath_valuations = "mmcsv/" + std::to_string(N_value) + "agents_valuations.csv";
        // Før du executer, lag mappe for dette runnet, så lag 3 csv.mapper inne i her igjen
        std::string filepath_executions = "TS/eqmaxi/executions/" + std::to_string(N_value) + "agents_executions.csv";
        std::string filepath_valuations = "TS/eqmaxi/valprofile/" + std::to_string(N_value) + "agents_valuations.csv";
        std::string filepath_solutions = "TS/eqmaxi/bestsolutions/" + std::to_string(N_value) + "best_solutions.csv";
        threads.emplace_back(run_algorithms, N_value, filepath_executions, filepath_valuations, filepath_solutions);
    }

    for (auto& thread : threads) {
        thread.join();
    }

    


    // ===================================================================
    // PARAMETERS for evolutionary algorithms
    /*
    const unsigned int T_SIZE = 8; // size for tournement selection
    const unsigned int VEC_SIZE = M; // Number of object variables in genotypes (number of elements in the vector)
    const unsigned int POP_SIZE = 25; // Size of population

    const unsigned int MAX_GEN = 250; // Maximum number of generation before STOP
    const unsigned int MIN_GEN = 100;  // Minimum number of generation before STOP
    const unsigned int STEADY_GEN = 70;  //stop after STEADY_GEN gen. without improvement

    const float P_CROSS = 0.7;      // Crossover probability
    const float P_MUT = 0.9;        // Mutation probability 

    const double EPSILON = 0.01;	// range for real uniform mutation
    double SIGMA = 0.3;	    	// std dev. for normal mutation
    
    const double hypercubeRate = 0.5;     // relative weight for hypercube Xover
    const double segmentRate = 0.5;  // relative weight for segment Xover

    const double uniformMutRate = 0.5;  // relative weight for uniform mutation
    const double detMutRate = 0.5;      // relative weight for det-uniform mutation
    const double normalMutRate = 0.5;   // relative weight for normal mutation


    // Create eval instances
    allocationEval<Allocation> fullEval;
    moFullEvalByCopy<moveNeighbor<Allocation>> moveEval(fullEval);

    // Evaluatig initial solution
    //fullEval(tmp1);



    // ############################
    // Initialisation of population
    // ############################
    // Based on a uniform generator
    
    eoUniformGenerator<size_t> uGen(0, N); // hadde N-1 her
    eoInitFixedLength<Allocation> random(VEC_SIZE, uGen);
    // Initialization of the population
    eoPop<Allocation> pop(POP_SIZE, random);

    // and evaluate it in one loop
    apply<Allocation>(fullEval, pop);       // STL syntax

    // OUTPUT
    pop.sort();
    
    std::cout << "Initial Population" << std::endl;
    std::cout << pop << std::endl;
    

    // ############################
    //          Selection
    // ############################
    //eoDetTournamentSelect<Allocation> selectOne(T_SIZE);

    

    eoProportionalSelect<Allocation> selectOne(pop);
    eoSelectPerc<Allocation> select(selectOne);

    eoGenerationalReplacement<Allocation> replace;
    //eoPlusReplacement<Allocation> replace;
    eoSegmentCrossover<Allocation> xoverS;

    eoHypercubeCrossover<Allocation> xoverA;

    eoPropCombinedQuadOp<Allocation> xover(xoverS, segmentRate);
    xover.add(xoverA, hypercubeRate);

    // ############################
    //          MUTATION
    // ############################

    
    // offspring(i) uniformly chosen in [parent(i)-epsilon, parent(i)+epsilon]
    eoUniformMutation<Allocation>  mutationU(EPSILON);
    
    // k (=1) coordinates of parents are uniformly modified
    eoDetUniformMutation<Allocation>  mutationD(EPSILON);
    // all coordinates of parents are normally modified (stDev SIGMA)
    //eoNormalMutation<Allocation>  mutationN(SIGMA);

    
    // Combine them with relative weights
    eoPropCombinedMonOp<Allocation> mutation(mutationU, uniformMutRate);
    mutation.add(mutationD, detMutRate);
    //mutation.add(mutationN, normalMutRate, true);

    // ############################
    //  CHECKPOINT
    // ############################
    // termination conditions: use more than one

    //STOP after MAX_GEN generations
    eoGenContinue<Allocation> genCont(MAX_GEN);
    // do MIN_GEN gen., then stop after STEADY_GEN gen. without improvement
    eoSteadyFitContinue<Allocation> steadyCont(MIN_GEN, STEADY_GEN);
    // stop when fitness reaches a target (here VEC_SIZE)
    //eoSteadyFitContinue<Allocation> fitCont(0);

    // do stop when one of the above says so
    eoCombinedContinue<Allocation> gaContinuator(genCont);
    gaContinuator.add(steadyCont);

    // The operators are  encapsulated into an eoTRansform object
    eoSGATransform<Allocation> transform(xover, P_CROSS, mutation, P_MUT);

    // ############################
    //          GENERATION
    // ############################
    //         the algorithm

    //eoEasyEA<Allocation> gga(gaContinuator, fullEval, select, transform, replace);
    // Apply algo on the population
    //gga(pop);

    
    //OUTPUT
    std::cout << "########################" << std::endl << std::endl;
    pop.sort();
    std::cout << "FINAL Population\n" << pop << std::endl;

    std::cout << "########################" << std::endl << std::endl;
    std::cout << "The best element: \n" << pop.best_element() << std::endl; */

}