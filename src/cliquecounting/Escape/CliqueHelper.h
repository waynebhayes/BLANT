#ifndef ESCAPE_CLIQUE_HELPER_H_
#define ESCAPE_CLIQUE_HELPER_H_

#include <algorithm>
#include <chrono>

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"
#include "Escape/Digraph.h"
#include "Escape/Utils.h"
#include "Escape/nCr.h"
#include "JointSort.h"

using namespace Escape;
using namespace std;
using namespace std::chrono;

struct VertexSet
{
    VertexIdx  nVertices;       //number of vertices in the VertexSet
    VertexIdx *vertices;        //array of vertices in the VertexSet

    bool contains(VertexIdx v) const;
    void print(FILE* f = stdout) const;
};

VertexSet* newVertexSet(VertexIdx nVertices)
{
    VertexSet *vs = new VertexSet();
    vs->nVertices = nVertices;
    vs->vertices = new VertexIdx[nVertices];
    return vs;
}

bool VertexSet::contains(VertexIdx v) const
{
    for (int i=0; i<nVertices; i++)
        if (vertices[i] == v)
            return true;
    return false;
}

void VertexSet::print(FILE* f) const
{
    if (nVertices == 0) fprintf(f, "Set empty\n");

    for (VertexIdx i = 0; i < nVertices; ++i)
        fprintf(f, "%lld   ", vertices[i]);
}

void delVertexSet(VertexSet* v)
{
    if (v != NULL)
    {
        if (v->nVertices != 0)
            delete[] v->vertices;
        delete v;
    }
}

VertexSet* intersect(VertexSet* v1, VertexSet* v2) 
{
    VertexSet *s=v2, *l=v1;
    if (v1->nVertices < v2->nVertices)
    {
        s = v1;
        l = v2;
    }

    VertexSet* v = newVertexSet(s->nVertices);
    int j = 0;
    for (int i=0; i<s->nVertices;i++)
    {
        if (l->contains(s->vertices[i]))
        {
            v->vertices[j] = s->vertices[i];
            j++;
        }
    }
    v->nVertices = j;
    return v;
}

struct PartialClique
{
    VertexSet* path;
    VertexSet* nbrs;
};

PartialClique* newPartialClique(VertexSet* path, VertexSet* nbrs)
{
    PartialClique* pc = new PartialClique;
    pc->path = path;
    pc->nbrs = nbrs;
    return pc;
}

void delPartialClique(PartialClique* pc)
{
    if (pc != NULL)
    {
        delVertexSet(pc->path);
        delVertexSet(pc->nbrs);
        delete pc;
    }
}

struct StackItem
{
    PartialClique *pc;
    StackItem *next;
};

void delStackItem(StackItem* si)
{
    if (si != NULL)
    {
        delPartialClique(si->pc);
        delete si;
    }
}


struct Stack
{
    StackItem *top;

    void push(PartialClique* pc);
    StackItem* pop();
    int i;
};

Stack* newStack()
{
    Stack* s = new Stack();
    s->top = NULL;
    s->i = 0;
    return s;
}

void delStack(Stack* s)
{
    if (s != NULL)
    {
        s->top = NULL;
        delete s;
    }
}

void Stack::push(PartialClique* pc)
{
    StackItem* si = new StackItem;
    si->pc = pc;
    si->next = top;
    top = si;
    i++;
}

StackItem* Stack::pop()
{
    StackItem *si = top;
    if (top == NULL) return NULL;
    top = si->next;
    si->next = NULL;
    i--;
    return si;
}

StackItem* newStackItem(PartialClique *pc, StackItem *next)
{
    StackItem* si = new StackItem;
    si->pc = pc;
    si->next = next;
    return si;
}

bool isClique(VertexIdx *vertices, int n, CGraph &CG)
{
    for (int i=0; i<n-1; i++)
    {
        for (int j=i+1; j<n; j++)
            if (CG.isEdge(vertices[i], vertices[j]) == -1)
                return false;
    }
    return true;
}

#endif
