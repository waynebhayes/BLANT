#include "blant-predict.h"
#include "blant-sampling.h"
#include "tinygraph.h"

#define MAX_NODES 28000 // may need to increase this for large networks

static int _predictOrd = -1;
int _predictOrbit1 = -1, _predictOrbit2 = -1;
static int _PredictMotifCount[MAX_K][MAX_K];
static float **_PredictCount; // will be allocated later
static GRAPH *_G;

// Count the specified orbits in the specified canonical motifs of g.
static void CountSubmotifs(TINY_GRAPH *g, int topCanon) {
  static int depth;
#if PARANOID_ASSERTS
  assert(g->n == _k && TinyGraphNumReachableNodes(g, 0) == g->n);
#endif
  int i, j, Gint = TinyGraph2Int(g, _k), GintOrd = _K[Gint];
  unsigned char perm[MAX_K];
  memset(perm, 0, _k);
  ExtractPerm(perm, Gint, false);
  if (depth == 0) {
    assert(topCanon == GintOrd); // we should only be called on canonicals
    for (i = 0; i < _k; i++)
      assert(perm[i] == i);
    for (i=0; i < _k; i++) for (j=0; j<_k; j++) _PredictMotifCount[i][j]=0;
  }
  // Check to see if the whole thing is the cononical of interest before we
  // start removing edges.
  if (GintOrd == _predictOrd) {
    for (i = 0; i < _k - 1; i++)
      for (j = i + 1; j < _k; j++)
        if (!TinyGraphAreConnected(g, perm[i], perm[j])) {
          if ((_orbitList[GintOrd][i] == _predictOrbit1 &&
               _orbitList[GintOrd][j] == _predictOrbit2) ||
              (_orbitList[GintOrd][i] == _predictOrbit2 &&
               _orbitList[GintOrd][j] == _predictOrbit1))
            ++_PredictMotifCount[(int)perm[i]][(int)perm[j]];
        }
  }

  // Now go about deleting edges recursively.
  for (i = 0; i < _k - 1; i++)
    for (j = i + 1; j < _k; j++) {
      if (TinyGraphAreConnected(g, i, j)) // if it's an edge, delete it.
      {
        TinyGraphDisconnect(g, i, j);
        if (TinyGraphNumReachableNodes(g, 0) == g->n) {
          ++depth;
          CountSubmotifs(g, topCanon);
          --depth;
        }
        TinyGraphConnect(g, i, j);
      }
    }
}

void Predict_ProcessGraphlet(GRAPH *G, unsigned Varray[], TINY_GRAPH *g,
                             Gint_type Gint, Gordinal_type GintOrd) {
#if PARANOID_ASSERTS
  assert(g->n == _k);
  assert(TinyGraphNumReachableNodes(g, 0) == g->n);
#endif
  unsigned i, j;
  assert(Gint == TinyGraph2Int(g, _k) && GintOrd == _K[Gint]);
  assert(_g_overcount);
  unsigned char perm[MAX_K];
  memset(perm, 0, _k);
  ExtractPerm(perm, Gint, false);
  TINY_GRAPH *gc = TinyGraphAlloc(g->n,false,false);
  Int2TinyGraph(gc, _canonList[GintOrd]);
  CountSubmotifs(gc, GintOrd);
  double totalWeight = 0.0;
  for (i = 0; i < _k; i++)
    totalWeight += 1.0 / GraphDegree(_G, Varray[i]);
  for (i = 0; i < _k - 1; i++)
    for (j = i + 1; j < _k; j++) {
      int g_u = perm[i], g_v = perm[j];
      int G_u = Varray[g_u], G_v = Varray[g_v];
      double remove_uv =
          1.0 / GraphDegree(_G, G_u) + 1.0 / GraphDegree(_G, G_v);
      _PredictCount[G_u][G_v] += (totalWeight - remove_uv)*(_PredictMotifCount[i][j] + _PredictMotifCount[j][i]);
    }
  TinyGraphFree(gc);
}

void Predict_Init(GRAPH *G) {
  assert(0 < _predictOrbit1 && _predictOrbit1 < _numOrbits);
  assert(0 < _predictOrbit2 && _predictOrbit2 < _numOrbits);

  if (_orbitCanonMapping[_predictOrbit1] != _orbitCanonMapping[_predictOrbit2])
    Fatal("orbit pair %d:%d must be in the same canonical, but they're in %d "
          "and %d, respectively",
          _predictOrbit1, _predictOrbit2, _orbitCanonMapping[_predictOrbit1],
          _orbitCanonMapping[_predictOrbit2]);

  _predictOrd = _orbitCanonMapping[_predictOrbit1];
  if (_predictOrbit2 < _predictOrbit1) {
    int tmp = _predictOrbit2;
    _predictOrbit2 = _predictOrbit1;
    _predictOrbit1 = tmp;
  }

  _G = G;
  _PredictCount = (float **)Calloc(G->n, sizeof(_PredictCount[0]));
  for (unsigned i = 0; i < G->n; i++)
    _PredictCount[i] = (float *)Calloc(G->n, sizeof(_PredictCount[i][0]));

  fprintf(
      stderr,
      "Sanity check: k %d, ordinal %d (integer %d), global orbitPair (%d,%d)\n",
      _k, _predictOrd, _canonList[_predictOrd], _predictOrbit1, _predictOrbit2);
}

void Predict_Flush(GRAPH *G) {
  // Currently a no-op; flushing is handled in Predict_Shutdown
}

void Predict_Shutdown(GRAPH *G) {
  unsigned i, j;
  fprintf(stderr, "count\tnode1\tnode2\n");
  for (i = 0; i < G->n - 1; i++)
    for (j = i + 1; j < G->n; j++) {
      float total = (_PredictCount[i][j] + _PredictCount[j][i]);
      if (total == 0)
        continue;
      printf("%g", total);
      if (_supportNodeNames)
        printf("\t%s\t%s\n", _nodeNames[i], _nodeNames[j]);
      else
        printf("\t%d\t%d\n", i, j);
    }
}
