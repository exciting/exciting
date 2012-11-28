template<class T>
class A {
public:
  static T x;
};

class B {
public:
  int y;
};

template <> B A<B>::x;

int main(void) {
  A<B>::x.y=1;
}
