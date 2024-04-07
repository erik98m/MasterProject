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

//-----------------------------------------------------------------------------
// Global problem description
// The number of agents, and the number of items
static int N, M;

// Valuations[n][m] is the utility of agent n and item m
static std::vector<std::vector<int>> Valuations;

//-----------------------------------------------------------------------------
// representation of solutions, and neighbors
// eoInt is based on eoVector class

// unsigned int for maximin and double with leximin
typedef eoInt<double> Allocation;

//-----------------------------------------------------------------------------
// fitness function and evaluation of solution
#include <eoEvalFunc.h>

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
        

        /*
        Maximin allocation

        // Find the agent with the lowest utility
        auto minUtility = std::numeric_limits<typename EOT::Fitness>::max();
        for (auto util : utility)
            minUtility = std::min(minUtility, util); */

        // Inform the allocation about its utility
        // fitness with leximin and minUtility with maximin
        allocation.fitness(fitness);
    }
};

//-----------------------------------------------------------------------------
// neighbor description

// Neighbor = How to compute the neighbor from the solution + information on it (i.e. fitness)
// all classes from paradisEO-mo use this template type

#include <neighborhood/moNeighbor.h>

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
#include <neighborhood/moNeighborhood.h>

// template <class EOT, class Fitness= unsigned int>
template <class EOT, class Fitness=typename EOT::Fitness>
class moveNeighborhood : public moNeighborhood<moveNeighbor<EOT, Fitness> >
{
public:
    typedef moveNeighbor<EOT, Fitness> Neighbor;

    /**
     * @return true if there is at least an available neighbor
     */
    virtual bool hasNeighbor(EOT& solution) {
        return N > 1 && M > 0;
    };

    /**
     * Initialization of the neighborhood
     * @param solution the solution to explore
     * @param current the first neighbor
     */
    virtual void init(EOT& solution, Neighbor& current) {
        // blir dette rikitg?
        nextItem = 0;
        nextAgent = 0;

        next(solution, current);
        //next(solution, current);
    }

    /**
     * Give the next neighbor
     * @param solution the solution to explore
     * @param current the next neighbor
     */
    virtual void next(EOT& solution, Neighbor& current) {
        current.setMove(nextItem, nextAgent);

        // Find out what the next neighbor should look like
        nextAgent++;

        // TODO: change to give random agent a random item?

        // Skip making a no-op move, and ensure correct wrapping
        while (nextItem < M) {
            if (nextAgent >= N) {
                nextAgent = 0;
                nextItem++;
            } else if (solution[nextItem] == nextAgent) {
                nextAgent++;
            } else break;
        }
    }

    /**
     * Test if there is again a neighbor
     * @param _solution the solution to explore
     * @return true if there is again a neighbor not explored
     */
    virtual bool cont(EOT& _solution) {
        return nextItem < M;
    }

    virtual std::string className() const {
        return "moveNeighborhood";
    }

private:
    // The pair of item and agent that makes the next neighbor
    int nextItem;
    int nextAgent;
};

//-----------------------------------------------------------------------------
// the simple Hill-Climbing local search
#include <algo/moSimpleHC.h>


int main()
{
    // TODO: Take in as parameters
    N = 5; // Number of agents
    M = 10; // Number of items

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0,N-1);

    Valuations.clear();
    Valuations.resize(N, std::vector<int>(M));
    // Assign random valuations to each item
    for (int n = 0; n < N; n++) {
        for (int m = 0; m < M; m++) {
            Valuations[n][m] = distr(gen);
        }
    }

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
    
    moveNeighbor<Allocation> move;
    moveNH.init(solution, move);
    moveEval(solution, move);
    move.print();
    while(moveNH.cont(solution)) {
        moveNH.next(solution, move);
        moveEval(solution, move);
        move.print();
    }

    // Define the simple Hill-Climbing 
    moSimpleHC<moveNeighbor<Allocation>> hc(moveNH, fullEval, moveEval);

    // apply the local search on the solution
    hc(solution);

    std::cout << "final: " << solution << std::endl ;

}