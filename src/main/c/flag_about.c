// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

/*! 
 * \file flag_about.c
 * Print information about the FLAG package, including version
 * and build numbers. 
 *
 * Usage: flag_about
 *
 */

#include <stdio.h>

int main(int argc, char *argv[]) {

  printf("%s\n", "==========================================================");
  printf("%s\n", "  FLAG package");
  printf("%s\n", "  3D Fourier-Laguerre transform on the Solid Sphere");
  printf("%s\n", "  By Boris Leistedt & Jason McEwen");

  printf("%s\n", "  See LICENSE.txt for license details.");

  printf("%s%s\n", "  Version: ", FLAG_VERSION);
  printf("%s%s\n", "  Build: ", FLAG_BUILD);
  printf("%s\n", "==========================================================");

  return 0;

}
