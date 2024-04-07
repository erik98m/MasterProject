#include <eval/moEval.h>
template< class Neighbor >
class moOneMaxIncrEval : public moEval<Neighbor>
{
public:

    typedef typename Neighbor::EOT EOT;

    /*
    * incremental evaluation of the neighbor for the oneMax problem
    * @param _solution the solution to move (bit string)
    * @param _neighbor the neighbor to consider (of type moBitNeigbor)
    */
    virtual void operator()(EOT & _solution, Neighbor & _neighbor) {
	if (_solution[_neighbor.index()] == 0)
	    _neighbor.fitness(_solution.fitness() + 1);
	else 
	    _neighbor.fitness(_solution.fitness() - 1);
    }
};
