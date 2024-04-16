#include "Cell.h"
#include "Disk.h"

Cell::Cell(int i_index)
{
    m_index = i_index;
    m_head_of_list = nullptr;
}

Cell::~Cell()
{

}

void Cell::set_head_of_list(Disk* i_head_of_list)
{
    m_head_of_list = i_head_of_list;
}

void Cell::add_neighbor(Cell* i_neighbor)
{
    m_neighbors.push_back(i_neighbor);
}

Disk* Cell::head_of_list()
{
    return m_head_of_list;
}

std::vector<Cell*> Cell::neighbors()
{
    return m_neighbors;
}
