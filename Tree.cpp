#include <map>
struct Node {
    int rank; // 0 for root
    int Parent_rank;
    int Sibling_rank;
    double *coefficients;
};
class Tree {
public:
    int N = 0; // Number of nodes in the tree
    int D; // Dimensionality of the tree
    void insert(Node node){ // define a new node
        Node * newnode = new Node;
        newnode->rank = node.rank;
        newnode->Parent_rank = node.Parent_rank;
        newnode->Sibling_rank = node.Sibling_rank;
        newnode->coefficients = new double[ncoefs];
        for (int i=0; i<ncoefs; i++) newnode->coefficients[i]=node.coefficients[i];
        nodemap[newnode->rank] = newnode;
        N++;
    };
    Node get(int rank) {// fetch a node
        Node node;
        node =
    };
    Tree(int Dims){
        D = Dims;
        ncoefs = 1;
        for (int i=0; i<Dims; i++) ncoefs*=2;
    };
    ~Tree(){
        for (int i=0; i<N; i++) delete[] (nodemap[i])->coefficients;
    };
private:
    std::map<int,Node*> nodemap;
    int ncoefs;
};
int main() {

}
