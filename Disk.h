#include "Eigen/Dense"

class Cell;

class Disk
{
    
private:
    int m_index;
    Cell* m_linkedCell;
    Disk* m_linkedDisk;
    
    double m_radius, m_mass, m_inertia;
    
    //translation
    Eigen::Vector3d m_r;
    Eigen::Vector3d m_v;
    Eigen::Vector3d m_a;
    Eigen::Vector3d m_F;

    //rotation
    double m_theta;
    Eigen::Vector3d m_w;
    Eigen::Vector3d m_alpha;
    Eigen::Vector3d m_M;

    //magnetic interaction
    Eigen::Vector3d m_mu;
    double m_chi; //angle between mu and e_x

    //material
    double E;
        
public:
    Disk(int,double,double,double,double,double,double,Eigen::Vector3d,double);
    ~Disk();
    
    void update_velocity(double);
    void update_position(double);
    void reset_force();
    void add_force(const Eigen::Vector3d&);
    void add_momentum(const Eigen::Vector3d&);
    void add_gravity_force(const Eigen::Vector3d&);
    void add_ddi_force(Eigen::Vector3d n, Disk& i_disk,Disk& j_disk);
    void set_linked_disk(Disk*);
    void set_linked_cell(Cell*);
    void print(int);
    
    Cell* linked_cell();
    Disk* linked_disk();

    //getteur
    int index();
    double radius();
    double mass();
    double theta();
    double get_E();
    double inertia();
    Eigen::Vector3d r();
    Eigen::Vector3d v();
    Eigen::Vector3d w();
    double x();
    double y();
    double vx();
    double vy();
    Eigen::Vector3d F();
    
    bool is_touching(double,double,double);

};
