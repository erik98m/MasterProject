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

// the general include for eo
#include <eo>
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

//-----------------------------------------------------------------------------
// Global problem description
// The number of agents, and the number of items
static int N, M;

// Valuations[n][m] is the utility of agent n and item m
static std::vector<std::vector<int>> Valuations;

//-----------------------------------------------------------------------------
// representation of solutions, and neighbors
// eoInt is based on eoVector class
// Change <unsigned int> for maximin and <double> with leximin
typedef eoInt<double> Allocation;

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

        // calculate fitnes with leximin
        
        // TODO: endre 0.1 til 1/(maxUti + 1)
        const double k = 0.1; 
        typename EOT::Fitness fitness = 0;

        // iterate from smallest to largest utility
        for (std::size_t i = 0; i < utility.size(); i++)
        {
            fitness += utility[i] * std::pow(k, i);
        }
        

        // Maximin allocation
    
        // Find the agent with the lowest utility
        auto minUtility = std::numeric_limits<typename EOT::Fitness>::max();
        for (auto util : utility)
            minUtility = std::min(minUtility, util);

        // Inform the allocation about its utility
        // fitness with leximin and minUtility with maximin
        allocation.fitness(fitness);
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

//-----------------------------------------------------------------------------

int main()
{
    // TODO: Take in as parameters
    N = 7; // Number of agents
    M = 8; // Number of items

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0,N-1);

    Valuations.clear();
    Valuations.resize(N, std::vector<int>(M));
    // Assign random valuations to each item
    /*
    for (int n = 0; n < N; n++) {
        for (int m = 0; m < M; m++) {
            Valuations[n][m] = distr(gen);
        }
    }*/

    // ####### Manually setup for testing purposes #######
    // Manually set valuations for the agents
    Valuations[0] = {1, 2, 3, 4, 5, 6, 7, 4}; // First agent
    Valuations[1] = {4, 3, 2, 1, 2, 5, 6, 7}; // Second agent
    Valuations[2] = {2, 1, 1, 4, 1, 1, 7, 7}; // Third agent
    Valuations[3] = {7, 7, 7, 7, 7, 0, 0, 1}; // Fourth agent
    Valuations[4] = {0, 0, 0, 6, 6, 6, 7, 7}; // Fifth agent
    Valuations[5] = {7, 6, 5, 4, 3, 2, 1, 0}; // Sixth agent
    Valuations[6] = {1, 1, 3, 1, 5, 6, 6, 7}; // Seventh agent
    
   /*
    Valuations[0] = {7, 7, 7, 7, 7, 7, 7, 7}; // First agent
    Valuations[1] = {7, 7, 7, 7, 7, 7, 7, 7}; // Second agents
    Valuations[2] = {7, 7, 7, 7, 7, 7, 7, 7}; // Third agents
    Valuations[3] = {7, 7, 7, 7, 7, 7, 7, 7}; // Fourth agents
    Valuations[4] = {7, 7, 7, 7, 7, 7, 7, 7}; // Fifth agents
    Valuations[5] = {7, 7, 7, 7, 7, 7, 7, 7}; // Sixth agents
    Valuations[6] = {7, 7, 7, 7, 7, 7, 7, 7}; // Seventh agents
    */
   
    // Initialize solution with random elements from 1 to n
    Allocation solution(M);
    for (int m = 0; m < M; m++)
    {
        solution[m] = distr(gen);
    }

    // Create eval instances
    allocationEval<Allocation> fullEval;
    moFullEvalByCopy<moveNeighbor<Allocation>> moveEval(fullEval);

    std::cout << "Initial Solution:" << std::endl;
    std::cout << solution << std::endl << std::endl;

    // move neighborhood
    moveNeighborhood<Allocation> moveNH;

    // Iterate through all neighbours and print their evaluation
    /*
    int moveCounter = 0;

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
     * the cooling schedule of the process
     *
     * ========================================================= */
    // Define the parameters for the cooling schedule
    double initialTemperature = 1;      // Initial temperature
    double coolingFactor = 0.9;        // Factor by which the temperature will be multiplied at each step
    unsigned span = 100;              // Number of iterations with the same temperature
    double finalTemperature = 0.01;  // Final temperature, stopping criterion

    moSimpleCoolingSchedule<Allocation> coolingSchedule(initialTemperature, coolingFactor, span, finalTemperature);

    /* =========================================================
     *
     * Checkpointing
     *
     * ========================================================= */
    moTrueContinuator<moveNeighbor<Allocation>> continuator;
    moCheckpoint<moveNeighbor<Allocation>> checkpoint(continuator);
    moCounterStat<Allocation> iterStat;
    checkpoint.add(iterStat);

    /* =========================================================
     *
     * Initialization and Execution of algorithms
     *
     * ========================================================= */
    //Define the simple Hill-Climbing 
    //moSimpleHC<moveNeighbor<Allocation>> hc(moveNH, fullEval, moveEval);
    // apply the local search on the solution
    //hc(solution);

    // Simulated Annealing
    moSA<moveNeighbor<Allocation>> localSearch1(moveNH, fullEval, moveEval, coolingSchedule, checkpoint);
    localSearch1(solution);

    // Printing out results
    std::cout << "final: " << solution << std::endl << std::endl ;
    std::cout << "Iterations: " << iterStat.getValue() << std::endl ;

}