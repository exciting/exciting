typedef Array<Real> PolytopeVertex;

typedef LinkedListIterator<PolytopeFace> LinkedListIteratorPolytopeFace;
typedef LinkedListIterator<PolytopeVertex> LinkedListIteratorPolytopeVertex;

class PolytopeFace {
 public:
  Array<Real> normal;
  Real c;
  Array<LinkedListIteratorPolytopeVertex> vertex;
  Array<LinkedListIteratorPolytopeFace> near_face;
};

