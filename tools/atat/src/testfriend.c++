class C {
public:
  friend int fr(int x) {
    return 1;
  }
};

int main(int argc, char *argv[]) {
  int y=fr(2);
}
