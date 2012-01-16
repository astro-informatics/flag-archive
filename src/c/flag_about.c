
#include <stdio.h>

int main(int argc, char *argv[]) {

  printf("%s\n", "==========================================================");

  printf("%s%s\n", "Version: ", SSHT_VERSION);
  printf("%s%s\n", "Build: ", SSHT_BUILD);
  printf("%s\n", "==========================================================");

  return 0;

}