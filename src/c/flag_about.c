
#include <stdio.h>

int main(int argc, char *argv[]) {

  printf("%s\n", "==========================================================");

  printf("%s%s\n", "Version: ", FLAG_VERSION);
  printf("%s%s\n", "Build: ", FLAG_BUILD);
  printf("%s\n", "==========================================================");

  return 0;

}