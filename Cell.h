#include <vector>

class Disk;

class Cell
{
    
private:
    int m_index;
    Disk* m_head_of_list;
    std::vector<Cell*> m_neighbors;
    
public:
    Cell(int);
    ~Cell();
    
    void set_head_of_list(Disk*);
    void add_neighbor(Cell*);
    
    Disk* head_of_list();
    std::vector<Cell*> neighbors();
};
