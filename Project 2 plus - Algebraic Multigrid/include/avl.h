#ifndef _AVL_H_
#define _AVL_H_

#include <vector>

struct avlNode{
    int val;
    int mxval;
    int size;
    avlNode *lson, *rson;
};

class avlTree{
private:
    avlNode *root;
    void build(avlNode* &rt, const int &l, const int &r, const std::vector<int> &vals);
    void increse(avlNode* rt, const int &k, const int &x);
    void setzero(avlNode* rt, const int &k);
    int getmaxID(avlNode* rt) const;

public:
    avlTree();
    void build(const std::vector<int> &vals);
    void increse(const int &k, const int &x);
    void setzero(const int &k);
    int getmaxID() const;
};

#endif